# -*- encoding: utf-8 -*-
'''
@File    :   parse.anchors.proali.py
@Time    :   2025/09/27 05:18:10
@Author  :   xiaodong li
@Version :   1.0.0
@Contact :   xiaodongli2405@gmail.com
'''

import sys

anchors = sys.argv[1]
out_block_info_file = sys.argv[2]

out_handle = open(out_block_info_file, 'w')
out_handle.write('BlockIndex\tRefChr\tRefStart\tRefEnd\tQueryChr\tQueryStart\tQueryEnd\tStrandt\tNumRecords\n')
with open(anchors) as f:
    next(f)
    next(f)

    block_index = 0
    ref_chr = ""
    ref_start = 0
    ref_end = 0

    query_chr = ""
    query_start = 0
    query_end = 0

    strand = ""

    flag = False
    for line in f:
        if line.startswith('#block begin'):
            block_index += 1
            num_records = 0
            ref_corr = set()
            query_corr = set()
            flag = True
            continue
        if line.startswith('#block end'):
            flag = False
            ref_start = min(ref_corr)
            ref_end = max(ref_corr)

            if strand == "+":
                query_start = min(query_corr)
                query_end = max(query_corr)
            else:
                query_start = max(query_corr)
                query_end = min(query_corr)

            out_handle.write(
                str(block_index) + '\t' + ref_chr + '\t' + str(ref_start) + '\t' + str(
                    ref_end) + '\t' + query_chr + '\t' + str(
                    query_start) + '\t' + str(query_end) + '\t' + strand + '\t' + str(num_records) + '\n')
            continue
        if flag:
            row_list = line.strip().split()
            ref_chr = row_list[0]
            query_chr = row_list[3]

            ref_start = int(row_list[1])
            ref_end = int(row_list[2])
            ref_corr.add(ref_start)
            ref_corr.add(ref_end)

            query_start = int(row_list[4])
            query_end = int(row_list[5])
            query_corr.add(query_start)
            query_corr.add(query_end)

            strand = row_list[6]

            num_records += 1


