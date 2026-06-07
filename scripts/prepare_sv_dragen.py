#!/usr/bin/env python3
"""
 * =========================================================
 * WGS/WES Germline Analysis Pipeline
 * =========================================================
 * Author   : Po-Yu Lin (林伯昱)
 * Institute: Department of Neurology and
 *            Department of Genomic Medicine,
 *            National Cheng Kung University Hospital
 * Contact  : p88124019@gs.ncku.edu.tw
 *
 * Copyright (c) 2026, Po-Yu Lin (林伯昱)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DISCLAIMER: This pipeline is provided "as is" without
 * warranty of any kind. The authors and their institution
 * make no representations or warranties regarding the
 * accuracy, completeness, or suitability of the analysis
 * results for any clinical or research purpose. Users are
 * solely responsible for validating and interpreting all
 * results.
 * =========================================================

prepare_sv_dragen.py
====================
目的：
    把 DRAGEN SV VCF（已過濾 PASS）的 INS 真實序列 ALT
    換成 <INS> symbolic allele。

    DRAGEN 的 INS ALT 可能長達數百至數千 bp，
    直接送 AnnotSV 會造成 TSV 欄位超長（Excel 無法顯示）。
    AnnotSV 對 symbolic allele <INS> 和真實序列的處理結果完全相同。

使用方式：
    python3 prepare_sv_dragen.py --input sv_pass.vcf.gz --output sv_filtered.vcf
"""

import argparse
import gzip
import sys


def convert_ins(input_path: str, output_path: str):
    n_ins_converted = 0
    n_other = 0

    opener = gzip.open if input_path.endswith('.gz') else open

    with opener(input_path, 'rt', encoding='utf-8') as fin, \
         open(output_path, 'w', encoding='utf-8') as fout:

        for line in fin:
            # header 行直接寫出
            if line.startswith('#'):
                fout.write(line)
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                fout.write(line)
                continue

            alt  = parts[4]
            info = parts[7]

            # 判斷是否為 INS 且 ALT 是真實序列（不是 symbolic allele）
            is_ins      = 'SVTYPE=INS' in info
            is_real_seq = not alt.startswith('<') and len(alt) > 10

            if is_ins and is_real_seq:
                parts[4] = '<INS>'
                fout.write('\t'.join(parts) + '\n')
                n_ins_converted += 1
            else:
                fout.write(line)
                n_other += 1

    print(f"[prepare_sv_dragen] INS 轉換：{n_ins_converted}，其他保留：{n_other}",
          file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description="把 DRAGEN SV VCF 的 INS 真實序列 ALT 換成 <INS>"
    )
    parser.add_argument('--input',  required=True, help='輸入 VCF（.vcf 或 .vcf.gz）')
    parser.add_argument('--output', required=True, help='輸出 VCF（未壓縮）')
    return parser.parse_args()


def main():
    args = parse_args()
    convert_ins(args.input, args.output)


if __name__ == '__main__':
    main()
