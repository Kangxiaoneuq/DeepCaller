import re
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import h5py
import pysam
import ctypes
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf

from collections import deque
import tensorflow_addons as tfa
from multiprocessing import Pool
from .utils import BedIndex, ref_encoding_dict


# Species → model folder name (must match main.py SPECIES_MAP)
SPECIES_MAP = {
    'potato':       '01_C88_Potato_default',
    'alfalfa':      '02_Bolivia_Alfalfa',
    'rose':         '03_Samantha_Rose',
    'sweetpotato':  '04_Tanzania_Sweetpotato_default',
    'syn_potato':   '05_SyntheticPotato_Potato',
}

# Mode → model subfolder prefix
MODE_PREFIX = {
    'speed':       '01_speed_model',
    'performance': '02_performance_mode',
}


def resolve_model_dir(species, mode):
    """
    Locate the model subdirectory for a given species and mode.
    Folder names follow the pattern: {prefix}_w={n}
    Returns (model_subdir_path, win_size).
    """
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
    project_root     = os.path.dirname(current_file_dir)
    species_dir      = os.path.join(project_root, "models", SPECIES_MAP[species])
    prefix           = MODE_PREFIX[mode]

    try:
        entries = os.listdir(species_dir)
    except FileNotFoundError:
        raise RuntimeError(f"Species model directory not found: {species_dir}")

    matched = [e for e in entries if e.startswith(prefix)]

    if not matched:
        raise RuntimeError(
            f"No model folder starting with '{prefix}' found in: {species_dir}"
        )

    folder = matched[0]
    match  = re.search(r'w=(\d+)', folder)

    if not match:
        raise RuntimeError(
            f"Cannot parse window size from folder name: {folder}"
        )

    win_size   = int(match.group(1))
    model_path = os.path.join(species_dir, folder)

    return model_path, win_size


def pileup_window_encoding(
    ref_window,
    seq_window,
    name_window,
    mq_window,
    strand_dict,
    win_size,
    bam_depth,
):
    """
    Encode a window of pileup information into a flat numerical array,
    split by strand direction and sequence features.
    """

    def encode_strand(seq_lst, names_lst, mq_lst, strand_dict, strand_type):

        strand_value = 1 if strand_type == 'pos' else 0

        mask = [strand_dict.get(name, True) == strand_value for name in names_lst]

        filtered_seqs = [seq for seq, m in zip(seq_lst, mask) if m]
        mq_lst_       = [mq  for mq,  m in zip(mq_lst,  mask) if m]

        snp_lst = ''.join(
            s[0].upper()
            for s in filtered_seqs
            if s != '*' and s and s[0] in 'ACGTacgt'
        )

        ins_lst = [s.split('+')[1].upper() for s in filtered_seqs if '+' in s]
        del_lst = [s.split('-')[1].upper() for s in filtered_seqs if '-' in s]

        base_count = [snp_lst.count(base) for base in 'ACGT']

        IS1 = pd.Series(ins_lst).value_counts().max() if ins_lst else 0
        IL1 = (
            int(re.match(r'^\d+', pd.Series(ins_lst).mode().iloc[0]).group(0))
            if ins_lst else 0
        )

        DS1 = pd.Series(del_lst).value_counts().max() if del_lst else 0
        DL1 = (
            int(re.match(r'^\d+', pd.Series(del_lst).mode().iloc[0]).group(0))
            if del_lst else 0
        )

        DR = filtered_seqs.count('*')

        result_mq = [[] for _ in range(7)]

        for idx, item in enumerate(filtered_seqs):
            if '*' in item:
                result_mq[6].append(idx)
                continue
            if not item or item[0] not in 'ACGTacgt':
                continue
            result_mq[['A', 'C', 'G', 'T'].index(item[0].upper())].append(idx)
            if IL1 and f"+{IL1}" in item:
                result_mq[4].append(idx)
            if DL1 and f"-{DL1}" in item:
                result_mq[5].append(idx)

        total_MQ = [sum(mq_lst_[idx] for idx in r) for r in result_mq]

        return base_count + [IS1, DS1, DR] + total_MQ

    results = []

    for i in range(1, win_size * 2 + 2):

        ref_encoding = [ref_encoding_dict.get(ref_window[i], 4)]

        pos_features = encode_strand(
            seq_window[i], name_window[i], mq_window[i], strand_dict, 'pos'
        )
        neg_features = encode_strand(
            seq_window[i], name_window[i], mq_window[i], strand_dict, 'neg'
        )

        results.append([bam_depth] + ref_encoding + pos_features + neg_features)

    return np.array(results, dtype=np.float32).flatten()


def select_encode_test(chrom, start, end, dct):

    fasta_path = dct['fasta_path']
    bam_path   = dct['bam_path']
    bed_path   = dct['bed_path']
    min_af    = dct['min_af']
    win_size   = dct['win_size']
    rd_floor   = dct['rd_floor']
    bam_depth  = dct['bam_depth']
    rd_divisor = dct['rd_divisor']
    max_id_len = dct['max_id_len']

    min_rd = max(round(int(bam_depth) / rd_divisor), rd_floor)

    fasta_file = pysam.FastaFile(fasta_path)
    chrom_len  = fasta_file.get_reference_length(chrom)
    bam_file   = pysam.AlignmentFile(bam_path, "rb", reference_filename=fasta_path)

    if bed_path:
        bed_idx = BedIndex(bed_path)

    strand_dict = {
        pread.qname: (pread.flag & 0x10) // 16
        for pread in bam_file.fetch(
            chrom, max(0, start - 100), end + 100, multiple_iterators=False
        )
    }

    ref_dict = {
        j: s.upper() if s in 'AGTC' else '*'
        for j, s in zip(
            range(max(1, start - 100), end + 101),
            fasta_file.fetch(chrom, max(1, start - 100) - 1, end + 100),
        )
    }

    pos_window, ref_window, seq_window, name_window, rd_window, mq_window = (
        deque(maxlen=win_size * 2 + 2) for _ in range(6)
    )

    enc_list, pos_list, ref_list, alt_list = [], [], [], []
    rd_list, ref_num_list, alt_num_list    = [], [], []

    for pcol in bam_file.pileup(
        chrom,
        max(1, start - win_size - 2),
        min(end + win_size, chrom_len),
        min_base_quality=0,
        truncate=True,
        multiple_iterators=False,
    ):
        pos_window.append(pcol.pos + 1)
        ref_window.append(ref_dict[pcol.pos + 1])
        seq_window.append(
            pcol.get_query_sequences(mark_matches=False, mark_ends=False, add_indels=True)
        )
        name_window.append(pcol.get_query_names())
        rd_window.append(pcol.get_num_aligned())
        mq_window.append(pcol.get_mapping_qualities())

        if len(pos_window) < win_size * 2 + 2:
            continue

        v_ref = ref_window[win_size + 1]
        v_rd  = rd_window[win_size + 1]

        if v_ref not in 'ACGT' or v_rd < min_rd:
            continue

        seq_lst = seq_window[win_size + 1]
        v_pos   = pos_window[win_size + 1]

        if bed_path and not bed_idx.contains(chrom, v_pos):
            continue

        snp_lst = [
            s[0].upper() for s in seq_lst if s != '*' and s[0].upper() != v_ref
        ]
        ins_lst = [s.split('+')[1].upper() for s in seq_lst if '+' in s]
        del_lst = [s.split('-')[1].upper() for s in seq_lst if '-' in s]

        snp_num = pd.Series(snp_lst).value_counts().max() if snp_lst else 0
        ins_num = pd.Series(ins_lst).value_counts().max() if ins_lst else 0
        del_num = pd.Series(del_lst).value_counts().max() if del_lst else 0

        alt_num  = max([snp_num, ins_num, del_num, 0])
        alt_freq = alt_num / v_rd

        if alt_freq < min_af:
            continue

        if snp_num >= ins_num and snp_num >= del_num:
            snp_element = pd.Series(snp_lst).mode().iloc[0] if snp_lst else None
            ref, alt = v_ref, snp_element

        elif ins_num >= snp_num and ins_num >= del_num:
            ins_element = (
                re.sub(r'^\d+', '', pd.Series(ins_lst).mode().iloc[0])
                if ins_lst else None
            )
            if len(ins_element) > max_id_len:
                continue
            ref, alt = v_ref, v_ref + ins_element

        else:
            del_length = (
                int(re.match(r'^\d+', pd.Series(del_lst).mode().iloc[0]).group(0))
                if del_lst else 0
            )
            if del_length > max_id_len:
                continue
            ref = ''.join(
                ref_dict[key]
                for key in range(v_pos, v_pos + del_length + 1)
                if key in ref_dict
            )
            alt = v_ref

        enc_list.append(
            pileup_window_encoding(
                ref_window, seq_window, name_window, mq_window,
                strand_dict, win_size, bam_depth,
            )
        )

        pos_list.append(np.uint32(v_pos))
        ref_list.append(ref)
        alt_list.append(alt)
        rd_list.append(np.uint16(v_rd))
        ref_num_list.append(np.uint16(sum(1 for s in seq_lst if s.upper() == v_ref)))
        alt_num_list.append(np.uint16(alt_num))

    df = pd.DataFrame({
        'chrom':   pd.Series([chrom] * len(pos_list), dtype='category'),
        'pos':     pd.Series(pos_list,     dtype='uint32'),
        'ref':     pd.Series(ref_list,     dtype='string'),
        'alt':     pd.Series(alt_list,     dtype='string'),
        'rd':      pd.Series(rd_list,      dtype='uint16'),
        'ref_num': pd.Series(ref_num_list, dtype='uint16'),
        'alt_num': pd.Series(alt_num_list, dtype='uint16'),
    })

    return enc_list, df


def process_chromosomes(dct):

    chrom         = dct['chrom']
    num_threads   = dct['num_threads']
    feature_dim   = dct['feature_dim']
    batch_size    = dct['batch_size']
    win_size      = dct['win_size']
    gpu_available = dct['gpu_available']
    model_path    = dct['model_path']
    stat_path     = dct['stat_path']
    work_dir      = dct['work_dir']

    seq_len = 2 * win_size + 1

    fasta_file = pysam.FastaFile(dct['fasta_path'])

    with h5py.File(stat_path, 'r') as f:
        mean = f['mean'][:]
        std  = f['std'][:]

    with tf.device('/CPU:0'):
        model = tf.keras.models.load_model(
            model_path,
            custom_objects={"F1Score": tfa.metrics.F1Score},
            compile=False,
        )

    print(f"[INFO] Encoding chromosome: {chrom}")

    chrom_len  = fasta_file.get_reference_length(chrom)
    chunk_size = chrom_len // num_threads

    chunks = [
        (chrom, start, min(start + chunk_size - 1, chrom_len), dct)
        for start in range(1, chrom_len + 1, chunk_size)
    ]

    chunks[-2] = (chrom, chunks[-2][1], chunks[-1][2], dct)
    chunks.pop()

    all_enc = []
    all_df  = pd.DataFrame()

    with Pool(processes=num_threads) as pool:
        for sub_enc, sub_df in pool.starmap(select_encode_test, chunks):
            all_enc.extend(sub_enc)
            sub_df = sub_df.dropna(axis=1, how='all')
            all_df = pd.concat([all_df, sub_df], axis=0, ignore_index=True)

    X_test  = np.array(all_enc)
    X_test -= mean
    X_test /= std
    X_test  = X_test.reshape(-1, seq_len, feature_dim)

    del all_enc
    ctypes.CDLL("libc.so.6").malloc_trim(0)

    device  = "/GPU:0" if gpu_available else "/CPU:0"
    y_preds = []

    with tf.device(device):
        print(f"[INFO] Running inference on {'GPU' if gpu_available else 'CPU'}...")
        for i in range(0, len(X_test), batch_size):
            y_preds.append(model.predict_on_batch(X_test[i:i + batch_size]))

    y_preds = np.concatenate(y_preds, axis=0)

    all_df["pred_prob"]  = np.max(y_preds,  axis=1).astype(np.float32)
    all_df["pred_label"] = np.argmax(y_preds, axis=1).astype(np.uint8)

    all_df.to_parquet(
        os.path.join(work_dir, f"{chrom}.parquet"),
        index=False,
        compression='zstd',
    )


def main():

    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('Required arguments')

    required.add_argument("--ref",           type=str,   required=True)
    required.add_argument("--bam",           type=str,   required=True)
    required.add_argument("--bed",           type=str,   default=None)
    required.add_argument("--chrom",         type=str,   required=True)
    required.add_argument("--ploidy",        type=int,   required=True)
    required.add_argument("--cpus",          type=int,   required=True)
    required.add_argument("--bam_depth",     type=int,   required=True)
    required.add_argument("--species",       type=str,   required=True)
    required.add_argument("--mode",          type=str,   required=True)
    required.add_argument("--min_af",        type=float, required=True)
    required.add_argument("--rd_floor",      type=int,   required=True)
    required.add_argument("--gpu_available",             required=True)
    required.add_argument("--work_dir",      type=str,   required=True)

    args = parser.parse_args()

    gpus = tf.config.list_physical_devices('GPU')

    if args.gpu_available.lower() == "true":
        if gpus:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
    else:
        tf.config.set_visible_devices([], 'GPU')

    model_path, win_size = resolve_model_dir(args.species, args.mode)
    stat_path            = os.path.join(model_path, "statistics.h5")
    model_path           = os.path.join(model_path, "model")

    dct = {
        'fasta_path'    : args.ref,
        'bam_path'      : args.bam,
        'bed_path'      : args.bed,
        'chrom'         : args.chrom,
        'ploidy'        : args.ploidy,
        'num_threads'   : args.cpus,
        'bam_depth'     : args.bam_depth,
        'win_size'      : win_size,
        'min_af'       : args.min_af,
        'rd_floor'      : args.rd_floor,
        'rd_divisor'    : args.ploidy,
        'gpu_available' : args.gpu_available.lower() == "true",
        'model_path'    : model_path,
        'stat_path'     : stat_path,
        'work_dir'      : args.work_dir,
        'batch_size'    : 8192,
        'feature_dim'   : 30,
        'max_id_len'    : 20,
    }

    process_chromosomes(dct)


if __name__ == "__main__":
    main()