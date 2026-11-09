#!/bin/python3.10
# -*- coding: utf-8 -*-
import sys
#sys.path.append("/drop/.local/lib/python3.10/site-packages")
from numba import njit, cuda, prange
import numpy as np
import pandas as pd

import math
import pyarrow.csv as pa_csv
from typing import Tuple
import argparse
import time
import os
import re
import gc

def read_meta_file(filepath: str) -> pd.DataFrame:

    parse_options = pa_csv.ParseOptions(delimiter="\t")
    read_options = pa_csv.ReadOptions(block_size=1e9)
    data = pa_csv.read_csv(
             filepath, parse_options=parse_options, read_options=read_options
    )
    data = data.to_pandas()
    return data


## CADD v1.6 RAW-score cut-offs equivalent to PHRED 15 / 20.
## These are the same values pre_pro2.sh / pre_1.sh use to build c5a (CADD15) and c5b (CADD20).
CADD_RAW_THRESHOLDS = {"15": 1.387112, "20": 2.097252}


def get_raw_scores(data: np.array) -> np.array:
    ## raw CADD columns, as read from the .meta file (blank -> NaN)
    return data[:, 16:25].astype("float")


def get_scores(data: np.array, min_raw: float = None) -> np.array:
    ##cadd score; reformat CADD range to 0-1
    ## the UKBB 200k cohort raw variant has CADD16 ranges from -18.793437 to 19.100986 including the prescore and permutated score (GRCh38_v1.6;VEP 100_GRCh38)
    ## `data` may be the full meta matrix or the raw-score block returned by get_raw_scores().
    ## min_raw: alleles with raw CADD below this are blanked (-> score 0), which is exactly
    ## what blanking them in c5a/c5b does upstream; this lets one CADD15 .meta give CADD20 too.
    scores = data[:, 16:25] if data.shape[1] > 9 else data
    scores = scores.astype("float").copy()
    if min_raw is not None:
        scores[scores < min_raw] = np.nan
    scores = (scores - (-19.811548)) / (25.028523 - (-19.811548))
    scores[np.isnan(scores)] = 0
    s0 = scores[:, 0]
    scores = np.insert(scores, 0, s0, axis=1)
    scores[np.isnan(scores)] = 0

    return scores

def get_allele_freq(data: np.array) -> np.array:
    ##allele frequency as it is in the UKBB cohort; this is currently based on the raw pVCF data.
    af_n = data[:, 6:15]
    af = np.where(
    (af_n == '') | (af_n == '0') | (af_n == '0.0') | (af_n.astype(float) == 0) | 
    (af_n.astype(str) == 'nan') | (af_n.astype(str) == 'NaN'),
    '3.86e-6',
    af_n)
    af = af.astype("float")
    ##the ref allele frequency
    af0 = 1 - np.nansum(af, axis=1)
    af0 = af0[:, np.newaxis]
    af[np.isnan(af)] = 1
    af = np.hstack([af0, af])

    return af

def read_karyo_samples(filepath: str) -> set:
    """
    Reads a list of sample IDs from a file (1 column, no header).
    These are samples with XO or XY karyotypes needing chrX genotype fixes.
    """
    with open(filepath) as f:
        samples = {line.strip() for line in f if line.strip()}
    return samples
def fix_chrX_gt(geno: str) -> str:
    """
    Fixes genotype strings like '0', '1', '.' (or '0:L', '1:P', '.:S') to diploid style. in chrX
    """
    if re.match(r"^0(?::.*)?$", geno):
        return "0/0"
    elif re.match(r"^1(?::.*)?$", geno):
        return "1/1"
    elif re.match(r"^2(?::.*)?$", geno):
        return "2/2"
    elif re.match(r"^3(?::.*)?$", geno):
        return "3/3"
    elif re.match(r"^4(?::.*)?$", geno):
        return "4/4"
    elif re.match(r"^5(?::.*)?$", geno):
        return "5/5"
    elif re.match(r"^6(?::.*)?$", geno):
        return "6/6"
    elif re.match(r"^7(?::.*)?$", geno):
        return "7/7"
    elif re.match(r"^8(?::.*)?$", geno):
        return "8/8"
    elif re.match(r"^9(?::.*)?$", geno):
        return "9/9"
    elif re.match(r"^10(?::.*)?$", geno):
        return "10/10"
    elif re.match(r"^\.(?::.*)?$", geno):
        return "./."
    return geno


def format_data(data: pd.DataFrame, xo_xy_samples: set) -> Tuple[np.array, np.array, np.array, np.array]:
    """
    Formats genotype matrix; applies chrX diploid transformation to samples with XO or XY karyotype.
    """
    # Extract metadata before converting to NumPy
    samples_header = data.columns[26:]
    chrom_column = data.iloc[:, 0].values
    samples_df = data.iloc[:, 26:]

    # Convert to NumPy
    data_np = data.to_numpy()
    samples_np = samples_df.to_numpy().astype(str)

    # Raw CADD block (thresholded per CADD cut-off later) and allele frequencies
    scores = get_raw_scores(data_np)
    af = get_allele_freq(data_np)

    # Normalize "0" → "0/0"
    samples_np[samples_np == "0"] = "0/0"

    # Get indices of samples present in xo_xy_samples
    fix_indices = [i for i, sample in enumerate(samples_header) if sample in xo_xy_samples]

    # Apply chrX fix only to rows with chromosome == "X"
    for row_idx, chrom in enumerate(chrom_column):
        if chrom == "chrX":
            for col_idx in fix_indices:
                samples_np[row_idx, col_idx] = fix_chrX_gt(samples_np[row_idx, col_idx])

    # Apply user-defined final fix
    fix_genotype_vec = np.vectorize(fix_genotype)
    samples_np = fix_genotype_vec(samples_np)
    samples_np = samples_np.astype("str")
    return scores, af, samples_np, samples_header

def fix_genotype(val):
    if re.fullmatch(r'(?:[\.0-9]+:[a-zA-Z]+|[01]|\.)', val):
        return '0/0'
    return val
    
def parse_command_line_args():
    parser = argparse.ArgumentParser(description="GenePy2 - Make score matrix")
    parser.add_argument("--gene", type=str, help="Gene name", required=True)
    parser.add_argument("--cadd", type=str, help="CADD score", required=True)
    parser.add_argument("--gpu", action="store_true", help="Use GPU")
    parser.add_argument(
        "--gpu-threads",
        type=int,
        help="Number of threads in each GPU block",
        nargs="?",
        const=256,
        default=256,
    )
    args = parser.parse_args()
    return args


@njit()
def get_score(S: np.array, af: np.array, db1: np.array, s_int: np.array):
    for i in range(db1.shape[0]):
        for j in prange(db1.shape[1]):
            if s_int[i, j, 0] == 0:
                db1[i][j] = S[i, s_int[i][j][1]] * (
                    -math.log10(af[i, s_int[i][j][0]] * af[i, s_int[i][j][1]])
                )

            elif s_int[i, j, 1] == 0:
                db1[i][j] = S[i, s_int[i][j][0]] * (
                    -math.log10(af[i, s_int[i][j][0]] * af[i, s_int[i][j][1]])
                )
            else:
                db1[i][j] = (
                    (S[i, s_int[i][j][0]] + S[i, s_int[i][j][1]])
                    * 0.5
                    * (-math.log10(af[i, s_int[i][j][0]] * af[i, s_int[i][j][1]]))
                )
    return db1


@cuda.jit
def get_score_kernel(S_d, db1_d, af_d, s_int_d):
    start = cuda.grid(1)
    stride = cuda.gridsize(1)

    for i in range(db1_d.shape[0]):
        for j in range(start, db1_d.shape[1], stride):

            if s_int_d[i, j, 0] == 0:
                db1_d[i][j] = S_d[i, s_int_d[i][j][1]] * (
                    -math.log10(af[i, s_int_d[i][j][0]] * af[i, s_int_d[i][j][1]])
                )

            elif s_int_d[i, j, 1] == 0:
                db1_d[i][j] = S_d[i, s_int_d[i][j][0]] * (
                    -math.log10(af[i, s_int_d[i][j][0]] * af[i, s_int_d[i][j][1]])
                )
            else:
                db1_d[i][j] = (
                    (S_d[i, s_int_d[i][j][0]] + S_d[i, s_int_d[i][j][1]])
                    * 0.5
                    * (-math.log10(af[i, s_int_d[i][j][0]] * af[i, s_int_d[i][j][1]]))
                )

def get_scores_cuda(
    scores: np.array,
    af: np.array,
    db1: np.array,
    samples_int: np.array,
    threads_per_block: int,
):

    num_cols = db1.shape[1]
    S_d = cuda.to_device(np.ascontiguousarray(scores))
    db1_d = cuda.to_device(np.ascontiguousarray(db1))
    af_d = cuda.to_device(np.ascontiguousarray(af))
    s_int_d = cuda.to_device(np.ascontiguousarray(samples_int))

    blockspergrid = (num_cols + threads_per_block - 1) // threads_per_block
    get_score_kernel[blockspergrid, threads_per_block](S_d, db1_d, af_d, s_int_d)

    S_d = None
    af_d = None
    s_int_d = None
    db1 = db1_d.copy_to_host()
    db1_d = None

    return db1


def index(data: np.array) -> np.array:
    int_data = data.ravel().view("uint8")[::4]
    int_data = int_data.reshape(data.shape[0], -1, 3) - 48
    int_data[(int_data[:, :, 2] > 253) | (int_data[:, :, 0] > 253)] = [0, 0, 0]
    int_data = int_data[:, :, [0, 2]].reshape(int_data.shape[0], int_data.shape[1], 2)
    return int_data.astype("uint8")
def nan_if(arr, value):
    return np.where(arr == value, np.nan, arr)
def score_db(
    samples: np.array,
    scores: np.array,
    af: np.array,
    gpu: bool,
    gpu_threads_per_block: int,
    samples_header=None,
    file_name: str = None,
    samples_int: np.array = None,
):
    ## samples_int can be passed in so the genotype index is built once and reused
    ## for every CADD cut-off.
    if samples_int is None:
        samples_int = index(samples)
    db1 = np.zeros_like(samples, dtype=float)
    if gpu:
        db1 = get_scores_cuda(
            scores, af, db1, samples_int, threads_per_block=gpu_threads_per_block
        )
    else:
        db1 = get_score(scores, af, db1, samples_int)
    db1[np.all(samples_int == np.asarray([0.0, 0.0], dtype="uint8"), axis=-1)] = 0.0
    out1 = np.nansum(nan_if(db1, "0.0"), axis=0)
    gg = np.array([file_name] * len(samples_header))
    U = np.vstack((samples_header, out1, gg)).T
    return U

def is_file_empty(file_path):
    return os.path.getsize(file_path) == 0
def is_file_empty_or_header_only(file_path):
    try:
        # Read only the first two lines
        df = pd.read_csv(file_path, sep='\t', nrows=2)
        
        # Check if there’s only a header or fewer than expected columns
        if df.shape[0] < 1:
            print(f"No data after header in file: {file_path}")
            return True
        if df.shape[1] < 26:
            print(f"Mismatch between header and expected columns in file: {file_path}")
            return True

    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return True
    print(f"correct file {file_path}")
    return False
    
def read_head(file_path, n=10):
    """
    Read the first n lines of a file.
    
    Args:
        file_path (str): Path to the text file.
        n (int): Number of lines to read (default is 10).
    
    Returns:
        list: List of the first n lines.
    """
    head_lines = []
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if i >= n:
                break
            head_lines.append(line.strip())  # strip() removes trailing newlines
    return head_lines
 
def output_meta_name(file_name: str, cadd: str) -> str:
    """ENSG..._CADD15.meta -> ENSG..._CADD<cadd>.meta (keeps the old per-threshold file names)."""
    m = re.search(r"_CADD\d+\.meta$", file_name)
    stem = file_name[: m.start()] if m else file_name[: -len(".meta")]
    return f"{stem}_CADD{cadd}.meta"


def parse_args():
    parser = argparse.ArgumentParser(
        description="GenePy: score one CADD15 .meta file at several CADD cut-offs in one pass"
    )
    parser.add_argument("meta", help="per-gene .meta file built at the lowest cut-off (CADD15)")
    parser.add_argument("kary", help="karyotype file (sample IDs with XO/XY karyotype)")
    parser.add_argument(
        "--cadd",
        nargs="+",
        default=["15", "20"],
        choices=sorted(CADD_RAW_THRESHOLDS),
        help="CADD PHRED cut-offs to output (default: 15 20)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    gene = args.meta
    file_name = os.path.basename(gene)
    xo_xy_samples = read_karyo_samples(args.kary)
    xo_xy_samples = ""  # NOTE: kept from the original script - this disables the chrX karyotype fix
    print(file_name)
    if not (file_name.startswith("ENSG") and file_name.endswith(".meta")):
        return

    if is_file_empty_or_header_only(gene):
        print(f"Skipping empty file: {gene}")
        ## (original script did not skip either; behaviour kept)

    # 1. read + parse genotypes ONCE (this is the expensive part)
    data = read_meta_file(filepath=gene)
    raw_scores, af, samples, samples_header = format_data(data, xo_xy_samples)
    del data
    gc.collect()
    samples_int = index(samples)

    # 2. score at each CADD cut-off, reusing the parsed genotypes
    for cadd in args.cadd:
        scores = get_scores(raw_scores, min_raw=CADD_RAW_THRESHOLDS[cadd])
        if (np.isnan(scores).sum()) < (scores.shape[0]):  # compute metascores if at least 1 variant
            out_name = output_meta_name(file_name, cadd)
            U1 = score_db(
                samples=samples,
                scores=scores,
                af=af,
                gpu=0,
                gpu_threads_per_block=0,
                samples_header=samples_header,
                file_name=out_name,
                samples_int=samples_int,
            )
            np.savetxt(out_name + ".txt", U1, fmt="%s", delimiter="\t")
            print(f"CADD{cadd}: wrote {out_name}.txt")

    del samples, samples_int, raw_scores, af, samples_header
    gc.collect()


if __name__ == "__main__":
    main()
