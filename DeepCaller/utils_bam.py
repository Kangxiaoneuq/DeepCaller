import os
import sys
import gzip
import pysam
import subprocess

def check_reference(ref_path):
    """
    Validate reference FASTA file and its index
    
    Args:
        ref_path: Path to reference FASTA file
        
    Raises:
        RuntimeError: If samtools is not available
        SystemExit: If reference file is missing or index generation fails
    """
   
    fai_path = f"{ref_path}.fai"
    if os.path.isfile(fai_path):
        print(f"[OK] Found FASTA index: {fai_path}")
        return
    
    print(f"[INFO] Generating FASTA index for {ref_path}...")
    try:
        subprocess.run(
            ["samtools", "faidx", ref_path],
            check=True,
            stderr=subprocess.PIPE
        )
        
        if not os.path.isfile(fai_path):
            raise RuntimeError("Index file not created")
        print(f"[OK] Generated FASTA index: {fai_path}")
        
    except subprocess.CalledProcessError as e:
        sys.exit(
            f"[ERROR] FASTA index generation failed:\n"
            f"{e.stderr.decode().strip()}"
        )

def get_chrom_list(ref_path, bed_path=None, user_input_chroms=None):
    """
    Get validated list of chromosomes to process.
    
    Args:
        ref_path: Path to reference FASTA file.
        bed_path: Path to optional BED file with specified chromosomes.
        user_input_chroms: List of user-specified chromosomes (None for all).
    
    Returns:
        List[str]: Validated chromosome names
    """
    # Step 1: Read reference FASTA file and get the list of chromosomes
    try:
        with pysam.FastaFile(ref_path) as fasta:
            ref_chroms = fasta.references
            if not ref_chroms:
                sys.exit(f"[ERROR] Empty reference: {ref_path}")
                
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read {ref_path}\nReason: {str(e)}")

    # Step 2: Handle the BED file if provided
    if bed_path:
        chrom_list = []
        seen = set()
        try:
            with (gzip.open(bed_path, "rt") if bed_path.endswith(".gz") else open(bed_path)) as f:
                for line in f:
                    if not line.strip() or line.startswith(("track", "browser", "#")):
                        continue
                    chrom = line.split("\t")[0]
                    if chrom not in seen:
                        chrom_list.append(chrom)
                        seen.add(chrom)
            
            if not chrom_list:
                sys.exit(f"[ERROR] No valid chromosomes found in BED file: {bed_path}")
        except Exception as e:
            sys.exit(f"[ERROR] Failed to read {bed_path}\nReason: {str(e)}")
        
        # Step 3: Ensure BED chromosomes are valid and in order
        invalid_bed_chroms = set(chrom_list) - set(ref_chroms)
        if invalid_bed_chroms:
            sys.exit(f"[ERROR] Invalid chromosomes in BED file: {invalid_bed_chroms}\nValid options: {' '.join(ref_chroms)}")

        print(f"* Use chromosomes from BED file: {chrom_list}")
        return chrom_list
    
    # Step 4: Handle user-specified chromosomes
    if user_input_chroms:
        invalid_chroms = set(user_input_chroms) - set(ref_chroms)
        if invalid_chroms:
            sys.exit(f"[ERROR] Invalid chromosomes:\n{invalid_chroms}\nValid options: {' '.join(ref_chroms)}")
        
        print(f"* Use user-specified chromosome(s): {user_input_chroms}")
        return user_input_chroms

    # Step 5: If no BED and no user input, return all reference chromosomes
    print("* No BED file or user-specified chromosomes, use all reference chromosomes")
    return list(ref_chroms)

                
def sort_bam(bam_path, num_threads, work_dir):
    """
    Verify BAM file is coordinate-sorted and properly indexed
    
    Args:
        bam_path: Path to input BAM file
        num_threads: Number of threads to use for indexing (if needed)
        work_dir: Directory for temporary files (unused in this version)
        
    Returns:
        str: Path to the input BAM file (if already sorted and indexed)
        
    Raises:
        SystemExit: If BAM is unsorted or indexing fails
    """

    # Check sort status
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            sort_order = bam.header.get("HD", {}).get("SO", "unspecified")
    except Exception as e:
        raise RuntimeError(f"Failed to read BAM file: {bam_path}\nReason: {str(e)}")

    # Reject unsorted BAMs
    if sort_order != "coordinate":
        sys.exit(
            f"[ERROR] BAM file is not coordinate-sorted: {bam_path}\n"
            "Please pre-sort the BAM with:\n"
            f"  samtools sort -@ {num_threads} -o sorted.bam {bam_path}\n"
            "Then re-run with the sorted file."
        )

    # Verify index
    bai_path = f"{bam_path}.bai"
    if not os.path.exists(bai_path):
        print(f"[INFO] Generating index for: {bam_path}...")
        try:
            subprocess.run(
                ["samtools", "index", "-@", str(num_threads), bam_path],
                check=True,
                stderr=subprocess.PIPE
            )
            print(f"[OK] Created index: {bai_path}")
        except subprocess.CalledProcessError as e:
            sys.exit(
                f"[ERROR] Index generation failed:\n"
                f"{e.stderr.decode().strip()}\n"
                "You can manually create it with:\n"
                f"  samtools index -@ {num_threads} {bam_path}"
            )
    else:
        print(f"[OK] Found BAM index: {bam_path}.bai")

    return bam_path

def filter_bam(args):
    """
    Filter BAM reads by quality criteria
    
    Args:
        args: Tuple containing (bam_path, chromosome, work_dir)
        
    Returns:
        str: Path to filtered BAM file
    """
    bam_path, chrom, work_dir = args
    
    temp_bam_path = os.path.join(work_dir, f"temp_{chrom}.bam")
    md_prefix     = os.path.join(work_dir, f"md_{chrom}")
    summary_path  = f"{md_prefix}.mosdepth.summary.txt"

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            with pysam.AlignmentFile(temp_bam_path, "wb", header=bam.header) as out_bam:
                for read in bam.fetch(chrom):
                    if (not read.is_unmapped and       #0x4
                        not read.mate_is_unmapped and  #0x8
                        not read.is_secondary and      #0x100
                        not read.is_qcfail and         #0x200
                        not read.is_duplicate and      #0x400
                        not read.is_supplementary and  #0x800
                        (not read.is_paired or read.is_proper_pair) and
                        read.mapping_quality >= 5):    
                        read.qname = f"{read.query_name}_{read.flag}"
                        out_bam.write(read)
                        
    except Exception as e:
        raise RuntimeError(f"Failed to filter {chrom} from {bam_path}\nReason: {str(e)}")
                            
    try:
        pysam.index(temp_bam_path)
    except Exception as e:
        raise RuntimeError(f"Indexing failed for {temp_bam_path}\nReason: {str(e)}")
    
    cmd = [
          "mosdepth",
          "-t", "1",
          "-x",  # fast-mode
          "-n",  # No per-base output 
          md_prefix,
          temp_bam_path
        ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"mosdepth failed for {chrom}\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDERR:\n{e.stderr}"
        ) 

def summarize_genome_depth(work_dir, chrom_list):
    total_cov = 0.0
    total_len = 0

    for chrom in chrom_list:
        prefix = os.path.join(work_dir, f"md_{chrom}")
        summary = f"{prefix}.mosdepth.summary.txt"
        if not os.path.exists(summary):
           continue

        length_i = None
        mean_i = None

        with open(summary, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("chrom"):
                    continue
                cols = line.split()
                if cols[0] == chrom or (cols[0].lower() == "total" and length_i is None):
                    try:
                        length_i = int(cols[1])
                        mean_i = float(cols[3])
                    except (IndexError, ValueError):
                        pass
                    
                    if cols[0] == chrom:
                        break

        if length_i and mean_i is not None:
            total_cov += mean_i * length_i
            total_len += length_i

    if total_len == 0:
        raise RuntimeError("No valid chromosome depth collected (total length is zero).")
 
    genome_mean = total_cov / total_len
    
    return round(genome_mean)