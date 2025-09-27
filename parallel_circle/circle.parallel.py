# -*- encoding: utf-8 -*-
'''
@File    :   circle.parallel.py
@Time    :   2025/09/27 05:18:10
@Author  :   xiaodong li
@Version :   1.0.0
@Contact :   xiaodongli2405@gmail.com
'''

import pandas as pd
import matplotlib.pyplot as plt
from math import pi
import numpy as np
import argparse


GAP_RATIO = 0.2
WIDTH = 10
HEIGHT = 10
DPI = 500
CHR_COLOR = "red"
INNER = 0.9

class Circle:
    def __init__(self, rf, qry, ars, fit, fig_bin):
        self.input_file = ars
        self.ref_length = rf
        self.query_length = qry
        self.flag = fit
        self.output_file_name = fig_bin

    @staticmethod
    def init_figure(width, height):
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif', 'Computer Modern Roman',
                                       'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                         'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter', 'serif']
        fig, ax = plt.subplots(figsize=(float(width), float(height)), facecolor='white')
        # fig.patch.set_alpha(0.5)
        ax.set_aspect('equal')
        return fig, ax

    @staticmethod
    def read_len(length_file, fit):
        _df = pd.read_csv(length_file, header=0, index_col=None, sep='\t')
        _df = _df[_df["chr"].str.endswith(fit)]
        return _df

    @staticmethod
    def get_circle_info(df_len):
        sp_total_chr_length = df_len['length'].sum()
        average_length = sp_total_chr_length / len(df_len)
        gap_length = int(average_length * GAP_RATIO)
        gap_nmb = len(df_len)-1
        circle_perimeter = sp_total_chr_length + gap_nmb * gap_length
        return circle_perimeter, gap_length

    @staticmethod
    def get_pos_list(df_dup, gap_length):
        chr_to_start = {}
        start_list = []
        end_list = []
        chr_length = df_dup['length'].tolist()
        chr_name = df_dup['chr'].tolist()
        i = 0
        for lgh in chr_length:
            if not start_list:
                start_list.append(0)
                end_list.append(lgh)
                chr_to_start[chr_name[i]] = 0
                i += 1
            else:
                start_list.append(end_list[i-1] + gap_length)
                end_list.append(start_list[i] + lgh)
                chr_to_start[chr_name[i]] = start_list[-1]
                i += 1
        chr_pos = list(zip(start_list, end_list))
        return chr_pos, chr_to_start

    @staticmethod
    def get_plt_chr_pos_x_y(chr_pos_pr_lt, sub_length_bp):
        divend = SUB_GENOME_RADIAN / sub_length_bp
        for pr in chr_pos_pr_lt:
            _x = []
            _y = []
            chr_st = pr[0]
            chr_end = pr[1]

            nrl_st = chr_st * divend
            nrl_end = chr_end * divend


            _x.append(np.sin(nrl_st))
            _y.append(np.cos(nrl_st))

            # get more dots
            length_list = np.linspace(nrl_st, nrl_end, 61800)
            for lgh in length_list:
                bt_x = np.sin(lgh)
                bt_y = np.cos(lgh)
                _x.append(bt_x)
                _y.append(bt_y)

            _x.append(np.sin(nrl_end))
            _y.append(np.cos(nrl_end))
            plt.plot(_x, _y, color=CHR_COLOR, zorder=3)
        return divend

    @staticmethod
    def get_query_plt_chr_pos_x_y(chr_pos_pr_lt, sub_length_bp):
        divend = SUB_GENOME_RADIAN / sub_length_bp
        for pr in chr_pos_pr_lt:
            _x = []
            _y = []
            chr_st = pr[0]
            chr_end = pr[1]

            nrl_st = chr_st * divend
            nrl_end = chr_end * divend


            _x.append(INNER * np.sin(nrl_st))
            _y.append(INNER * np.cos(nrl_st))

            # get more dots
            length_list = np.linspace(nrl_st, nrl_end, 61800)
            for lgh in length_list:
                bt_x = np.sin(lgh)
                bt_y = np.cos(lgh)
                _x.append(INNER * bt_x)
                _y.append(INNER * bt_y)

            _x.append(INNER * np.sin(nrl_end))
            _y.append(INNER * np.cos(nrl_end))
            plt.plot(_x, _y, color=CHR_COLOR, zorder=3)
        return divend

    @staticmethod
    def bezier3(a, b, c, d, t):
        fvalue = []
        for i in t:
            y = a * (1 - i) ** 3 + 3 * b * i * (1 - i) ** 2 + 3 * c * (1 - i) * i ** 2 + d * i ** 3
            fvalue.append(y)
        return fvalue

    def plot_closed_region(self, pos1_radian, pos2_radian, pos3_radian, pos4_radian, radius1, radius2, clr, zorder):
        ratio = 0.5

        t = np.arange(pos1_radian, pos2_radian, pi / 61800)
        y = list(radius1 * np.cos(t))
        x = list(radius1 * np.sin(t))

        pos2_x = radius1 * np.sin(pos2_radian)
        pos2_y = radius1 * np.cos(pos2_radian)
        pos4_x = radius2 * np.sin(pos4_radian)
        pos4_y = radius2 * np.cos(pos4_radian)
        t = np.arange(0, 1.01, 0.01)
        dx = pos4_x - pos2_x
        dy = pos4_y - pos2_y
        p1, p2, p3, p4 = pos2_x, dx * ratio + pos2_x, -dx * ratio + pos4_x, pos4_x
        x1 = self.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos2_y, (1 - ratio) * dy + pos2_y, -(1 - ratio) * dy + pos4_y, pos4_y
        y1 = self.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        t = np.arange(pos4_radian, pos3_radian, pi / 61800)
        x1 = list(radius2 * np.sin(t))
        y1 = list(radius2 * np.cos(t))
        x = x + x1
        y = y + y1

        pos3_x = radius2 * np.sin(pos3_radian)
        pos3_y = radius2 * np.cos(pos3_radian)
        pos1_x = radius1 * np.sin(pos1_radian)
        pos1_y = radius1 * np.cos(pos1_radian)
        t = np.arange(0, 1.01, 0.01)

        dx = pos1_x - pos3_x
        dy = pos1_y - pos3_y
        p1, p2, p3, p4 = pos3_x, dx * ratio + pos3_x, -dx * ratio + pos1_x, pos1_x
        x1 = self.bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos3_y, (1 - ratio) * dy + pos3_y, -(1-ratio) * dy + pos1_y, pos1_y
        y1 = self.bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        x.append(x[0])
        y.append(y[0])
        plt.fill(x, y, alpha=0.75, color=clr, zorder=zorder, linewidth=0)
        return x, y

    def plot_block(self, anchors, rf_mp, qry_mp, r_r, q_r, flg):
        with open(anchors, 'r') as f:
            next(f)
            block = []
            for line in f:
                line_list = line.strip().split()

                ref_chr = line_list[1]
                query_chr = line_list[4]

                if not ref_chr.endswith(flg) or not query_chr.endswith(flg):
                    continue
                r_add_value = rf_mp[ref_chr]
                ref_start = int(line_list[2]) + r_add_value
                ref_end = int(line_list[3]) + r_add_value


                q_add_value = qry_mp[query_chr]
                query_start = int(line_list[5]) + q_add_value
                query_end = int(line_list[6]) + q_add_value


                direction = line_list[7]

                block.append([])
                block[-1].append(ref_start * r_r)
                block[-1].append(ref_end * r_r)
                block[-1].append(query_start * q_r)
                block[-1].append(query_end * q_r)
                block[-1].append(direction)

        for blk in block:
            if blk[4] == "+":
                color = "#F0F0F0"
                zorder = 1
            else:
                color = "#66AD56"
                zorder = 2
            self.plot_closed_region(float(blk[0]), float(blk[1]), float(blk[2]), float(blk[3]), 1, INNER, color, zorder)

    def run(self):
        fig, ax = self.init_figure(WIDTH, HEIGHT)

        ref_df = self.read_len(self.ref_length, self.flag)
        query_df = self.read_len(self.query_length, self.flag)

        # get reference all chromosome's total length
        ref_circle_perimeter, ref_gap_length = self.get_circle_info(ref_df)
        # plot ref chr
        ref_chr_pos, ref_chr_to_start = self.get_pos_list(ref_df, ref_gap_length)
        r_div = self.get_plt_chr_pos_x_y(ref_chr_pos, ref_circle_perimeter)

        # get query all chromosome's total length
        query_circle_perimeter, query_gap_length = self.get_circle_info(query_df)
        # plot query chr
        query_chr_pos, query_chr_to_start = self.get_pos_list(query_df, ref_gap_length)
        q_div = self.get_query_plt_chr_pos_x_y(query_chr_pos, query_circle_perimeter)

        _ = self.plot_block(anchors, ref_chr_to_start, query_chr_to_start, r_div, q_div, self.flag)

        plt.axis('off')
        plt.savefig(self.output_file_name, dpi=DPI, bbox_inches='tight', transparent=False)
        return "I love you!"



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='anchors file from AnchorWave proali subcommand.',
        prog="parallel circle")
    subparsers = parser.add_subparsers(title='Anchors file visualization from AnchorWave proali subcommand', dest='analysis')

    parser_sub_line_proali = subparsers.add_parser('parallel_plot',
                                                   help='Anchors file from AnchorWave proali to visualization.')

    parser_sub_line_proali.add_argument('-a', '--anchors', dest='anchors',
                                        help="anchors file.", metavar="")
    parser_sub_line_proali.add_argument('-r', '--ref', dest='ref', help="ref length file.",
                                        metavar="")
    parser_sub_line_proali.add_argument('-q', '--query', dest='query',
                                        help="query length file.",
                                        metavar="")
    parser_sub_line_proali.add_argument('-f', '--flag', dest='flag',
                                        help="sub_genome (A or B or D)", metavar="")
    parser_sub_line_proali.add_argument('-o', '--out_file', dest='out_file',
                                        help="Figure name", metavar="")
    parser_sub_line_proali.add_argument('-angle', '--angle', dest='angle',
                                        help="angle(0-360)", metavar="")

    args = parser.parse_args()

    # chinese spring wheat
    ref = args.ref
    # query jm47
    query = args.query
    # anchors
    anchors = args.anchors
    # flag (A or B or D)
    flag = args.flag
    # out
    out_file = args.out_file
    #radian
    SUB_GENOME_RADIAN = float(args.angle) / 360 * 2 * pi
    test_class = Circle(ref, query, anchors, flag, out_file).run()
    print(test_class)
