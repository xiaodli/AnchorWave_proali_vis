# -*- encoding: utf-8 -*-
'''
@File    :   line_proali_pangenome.py
@Time    :   2025/03/23 21:22:45
@Author  :   xiaodong li
@Version :   1.0.0
@Contact :   xiaodongli2405@gmail.com
'''


import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.transforms import Bbox
import pandas as pd
from collections import OrderedDict
import numpy as np
from math import pi
import sys
import os
# import time


class FileEmptyError(Exception):
    pass

def file_empty(file_path):
    try:
        file_path = os.path.abspath(file_path)
        if os.path.exists(file_path):
            if os.path.getsize(file_path) > 0:
                pass
            else:
                raise FileEmptyError
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        exist_error_message = f"{file_path} don't exist or the path physically exists but permission is not granted to execute os.stat() on the requested file."
        print(exist_error_message)
        sys.exit(1)
    except FileEmptyError:
        empty_error_message = "{0} is empty.".format(file_path)
        print(empty_error_message)
        sys.exit(1)
    except OSError:
        ose_error_message = f"{file_path} does not exist or is inaccessible."
        print(ose_error_message)
        sys.exit(1)

def bezier3(a, b, c, d, t):
    fvalue = []
    for i in t:
        y = a * (1-i)**3 + 3 * b * i * (1-i)**2 + 3 * c * (1-i) * i**2 + d * i**3
        fvalue.append(y)
    return fvalue


class Line:
    def __init__(self, config_pra, parameter):
        # gaps between chromosome, chr:gap = 4: 1
        self.gap_ratio = 6
        self.ref_height = 0.3
        self.height_gap = 0.07
        self.query_height = 0.37
        self.width = 0.0015
        self.dpi = 2500
        self.sp_chr_color_comma_sep = ""

        self.chr_font_size = 7
        self.species_name_font_size = 7
        self.hide_chr = False
        self.italic = False
        self.figsize = "14,14"
        self.color_style = "four_colors"
        self.actual_len = False
        self.gap_style = "compact"

        for i in config_pra.sections():
            if i == 'line':
                for key in config_pra[i]:
                    setattr(self, key, config_pra[i][key])

        for key, value in vars(parameter).items():
            if key != "func" and key != "analysis" and value is not None:
                setattr(self, key, value)

    def determine_fig_par(self, len_number):
        if len_number == 2:
            x = 0.4 / (6 * len_number - 5)
        else:
            x = 0.7 / (6 * len_number - 5)
        self.height_gap = 5 * x
        self.query_height = 0.3 + 5 * x

    def get_sp_length(self):
        length_file_list = self.length_file.split(",")
        sp_length = []
        i = 0
        for le in length_file_list:
            le = le.strip()
            if len(le) == 0:
                continue
            file_empty(le)
            length = pd.read_csv(le, sep='\t', header=0, index_col=None)
            length['chr'] = length['chr'].astype(str)
            length['length'] = length['length'].astype(int)
            sp_length.append(length['length'].sum())
            i += 1
        return sp_length

    @staticmethod
    def get_expand_ratio(length_list):
        max_length = max(length_list)
        ratio_list = []
        for i in length_list:
            if i < max_length * 0.5:
                ratio_list.append(max_length * 0.5 / i)
            elif max_length * 0.5 <= i < max_length * 0.65:
                ratio_list.append(max_length * 0.65 / i)
            elif max_length * 0.65 <= i < max_length * 0.8:
                ratio_list.append(max_length * 0.8 / i)
            else:
                ratio_list.append(1)
        return ratio_list

    @staticmethod
    def read_len(file, prefix):
        length = pd.read_csv(file, sep='\t', header=0, index_col=None)
        length['chr'] = length['chr'].astype(str)
        length['chr'] = prefix + length['chr']
        length['length'] = length['length'].astype(int)
        return length

    def split_length_adjust(self, strip_prefix, ratio_list):
        length_file_list = self.length_file.split(",")
        rm_blank = []
        for idx, le in enumerate(length_file_list):
            le = le.strip()
            if len(le) == 0:
                continue
            rm_blank.append(le)
        strip_length_file_list_df = []

        for idx, le in enumerate(rm_blank):
            file_empty(le)
            ratio = ratio_list[idx]
            length = self.read_len(le, strip_prefix[idx])
            length['length'] = length['length'] * ratio
            strip_length_file_list_df.append(length)
        return strip_length_file_list_df

    @staticmethod
    def split_conf(conf, separator):
        new_conf = []
        split_lt = conf.split(separator)
        for ele in split_lt:
            ele = ele.strip()
            if len(ele) == 0:
                continue
            new_conf.append(ele)
        return new_conf

    def get_sp_chr_color_dict(self, sp_list):
        if self.sp_chr_color_comma_sep:
            color_list = self.split_conf(self.sp_chr_color_comma_sep, ",")
            i = 0
            while len(color_list) < len(sp_list):
                dup_color = color_list[i]
                color_list.append(dup_color)
                if i >= len(color_list)-1:
                    i = 0
                else:
                    i += 1
        else:
            color_list = ["black"] * len(sp_list)
        sp_chr_color_dict = dict(list(zip(sp_list, color_list)))
        return sp_chr_color_dict

    def get_my_cmap(self):
        if self.color_style == "rainbow":
            cmap = self.set_rainbow()
        elif self.color_style == "husl":
            cmap = self.set_husl_palette()
        elif self.color_style == "four_colors":
            cmap = self.set_palette()
        else:
            cmap = self.set_husl_palette()
        return cmap

    @staticmethod
    def set_palette():
        colors = ['#FFA07A', '#4b6870', '#4169E1', '#9370DB']
        cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
        return cmap

    @staticmethod
    def set_husl_palette():
        colors = sns.husl_palette(as_cmap=True)
        cmap = LinearSegmentedColormap.from_list("husl_256", colors(np.linspace(0, 1, 256)))
        return cmap

    @staticmethod
    def set_rainbow():
        colors = plt.get_cmap("gist_rainbow")
        cmap = LinearSegmentedColormap.from_list("rainbow", colors(np.linspace(0, 1, 256)))
        return cmap

    @staticmethod
    def circle(radius, x):
        if x > 0.8:
            print(x)
        y = np.sqrt(radius ** 2 - x ** 2)
        return y

    @staticmethod
    def line_get_pos_list_compact(df_dup, gap_length, total_length, flag):
        tmp_start = (total_length - df_dup['length'].sum() - len(df_dup) * gap_length) / 2
        chr_to_start = {}
        start_list = []
        end_list = []
        chr_length = df_dup['length'].tolist()
        chr_name = df_dup['chr'].tolist()
        i = 0
        for lgh in chr_length:
            if not start_list:
                if flag:
                    start_list.append(gap_length + tmp_start)
                    end_list.append(gap_length + tmp_start + lgh)
                    chr_to_start[chr_name[i]] = gap_length + tmp_start
                else:
                    start_list.append(gap_length)
                    end_list.append(gap_length + lgh)
                    chr_to_start[chr_name[i]] = gap_length
                i += 1
            else:
                start_list.append(end_list[i-1] + gap_length)
                end_list.append(start_list[i] + lgh)
                chr_to_start[chr_name[i]] = start_list[-1]
                i += 1
        return start_list, end_list, chr_to_start

    @staticmethod
    def line_get_pos_list(df_dup, gap_length):
        chr_to_start = {}
        start_list = []
        end_list = []
        chr_length = df_dup['length'].tolist()
        chr_name = df_dup['chr'].tolist()
        i = 0
        for lgh in chr_length:
            if not start_list:
                start_list.append(gap_length)
                end_list.append(gap_length + lgh)
                chr_to_start[chr_name[i]] = gap_length
                i += 1
            else:
                start_list.append(end_list[i-1] + gap_length)
                end_list.append(start_list[i] + lgh)
                chr_to_start[chr_name[i]] = start_list[-1]
                i += 1
        return start_list, end_list, chr_to_start

    @staticmethod
    def read_proali_collinearity(qry_prefix, ref_prefix, collinearity, chr_list, chr_to_start):
        # left -> ref;  right -> query
        data = []
        ref_chr_list = []
        query_chr_list = []
        direction_list = []
        block_index = 0
        block = []
        with open(collinearity) as f:
            print("read", collinearity, "....")
            _ = next(f)
            _ = next(f)
            flag = True
            flag_number = True
            for line in f:
                if line.startswith("#block begin"):
                    flag = True
                    flag_number = True
                    if block:
                        data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
                        block = []
                elif line.startswith("#block end"):
                    continue
                else:
                    if flag:
                        chr_pair = [line.split()[0], line.split()[3]]
                        if ref_prefix + chr_pair[0] not in chr_list or qry_prefix + chr_pair[1] not in chr_list:
                            flag = False
                        else:
                            if flag_number:
                                flag = True
                                ref_chr = ref_prefix + chr_pair[0]
                                ref_chr_list.append(ref_chr)
                                query_chr = qry_prefix + chr_pair[1]
                                query_chr_list.append(query_chr)
                                direction_list.append(line.split()[6])
                                block_index += 1
                                flag_number = False
                        if flag:
                            line_list = line.split()
                            block.append([chr_to_start[ref_chr] + int(line_list[2]), chr_to_start[query_chr] + int(line_list[5])])
                    else:
                        continue
            data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
            print("parse", collinearity, "success")
        return data, ref_chr_list, query_chr_list, direction_list

    @staticmethod
    def read_proali_collinearity_ratio(qry_prefix, ref_prefix, collinearity, chr_list, chr_to_start, ref_query_ratio):
        ref_ratio = ref_query_ratio[0]
        query_ratio = ref_query_ratio[1]
        # left -> ref;  right -> query
        data = []
        ref_chr_list = []
        query_chr_list = []
        direction_list = []
        block_index = 0
        block = []
        with open(collinearity) as f:
            print("read", collinearity, "....")
            _ = next(f)
            _ = next(f)
            flag = True
            flag_number = True
            for line in f:
                if line.startswith("#block begin"):
                    flag = True
                    flag_number = True
                    if block:
                        data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
                        block = []
                elif line.startswith("#block end"):
                    continue
                else:
                    if flag:
                        chr_pair = [line.split()[0], line.split()[3]]
                        if ref_prefix + chr_pair[0] not in chr_list or qry_prefix + chr_pair[1] not in chr_list:
                            flag = False
                        else:
                            if flag_number:
                                flag = True
                                ref_chr = ref_prefix + chr_pair[0]
                                ref_chr_list.append(ref_chr)
                                query_chr = qry_prefix + chr_pair[1]
                                query_chr_list.append(query_chr)
                                direction_list.append(line.split()[6])
                                block_index += 1
                                flag_number = False
                        if flag:
                            line_list = line.split()
                            block.append([chr_to_start[ref_chr] + int(line_list[2]) * ref_ratio, chr_to_start[query_chr] + int(line_list[5]) * query_ratio])
                    else:
                        continue
            data.append([block[0][0], block[0][1], block[-1][0], block[-1][1]])
            print("parse", collinearity, "success")
        return data, ref_chr_list, query_chr_list, direction_list

    @staticmethod
    def set_husl_palette():
        colors = sns.husl_palette(as_cmap=True)
        cmap = LinearSegmentedColormap.from_list("husl_256", colors(np.linspace(0, 1, 256)))
        return cmap

    @staticmethod
    def set_palette():
        colors = ['#FFA07A', '#4b6870', '#4169E1', '#9370DB']
        cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)
        return cmap

    @staticmethod
    def plot_line_chr(start, end, width, height):
        new_start = start + width / 2
        new_end = end - width / 2

        # right
        t = np.arange(-pi/2, pi/2, pi / 180)
        x = list(width / 2 * np.cos(t) + new_end)
        y = list(width / 2 * np.sin(t) + height + width / 2)

        # top
        x.append(new_end)
        y.append(height + width)
        x.append(new_start)
        y.append(height + width)

        # left
        t = np.arange(pi/2, 3*pi/2, pi / 180)
        x1 = list(width / 2 * np.cos(t) + new_start)
        y1 = list(width / 2 * np.sin(t) + height + width / 2)
        x += x1
        y += y1

        # bottom
        x.append(new_start)
        y.append(height)
        x.append(new_end)
        y.append(height)

        return x, y

    def plot_collinearity_region(self, pos1, pos2, pos3, pos4, query_height, ref_height, jg_qr_st, jg_qr_ed, jg_rf_st, jg_rf_ed):
        ratio = 0.316
        # ref_mar = ref_height + self.width * 1.1
        # query_mar = query_height - self.width * 0.1
        ref_mar = ref_height + self.width * 1.0
        query_mar = query_height
        x, y = [], []
        # p1->p3(ref_block)
        t = np.linspace(pos1, pos3, 1000)
        for i in t:
            if jg_rf_st - self.width / 2 <= i < jg_rf_st:
                will_x = i
                will_y = ref_height + self.width / 2 + self.circle(self.width / 2, jg_rf_st - i)
                x.append(will_x)
                y.append(will_y)
            elif jg_rf_st <= i <= jg_rf_ed:
                will_x = i
                will_y = ref_mar
                x.append(will_x)
                y.append(will_y)
            else:
                will_x = i
                will_y = ref_height + self.width / 2 + self.circle(self.width / 2, i - jg_rf_ed)
                x.append(will_x)
                y.append(will_y)

        # p3->p4
        pos3_x = x[-1]
        pos3_y = y[-1]
        if jg_qr_st - self.width / 2 < pos4 < jg_qr_st:
            will_x = pos4
            will_y = query_height + self.width / 2 - self.circle(self.width / 2, jg_qr_st - pos4)
        elif jg_qr_st <= pos4 <= jg_qr_ed:
            will_x = pos4
            will_y = query_mar
        else:
            will_x = pos4
            will_y = query_height + self.width / 2 - self.circle(self.width / 2, pos4 - jg_qr_ed)
        t = np.arange(0, 1.01, 0.005)
        dx = will_x - pos3_x
        dy = will_y - pos3_y
        p1, p2, p3, p4 = pos3_x, dx * ratio + pos3_x, -dx * ratio + will_x, will_x
        x1 = bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos3_y, (1 - ratio) * dy + pos3_y, -(1 - ratio) * dy + will_y, will_y
        y1 = bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        # p4->p2
        t = np.linspace(pos4, pos2, 1000)
        for i in t:
            if jg_qr_st - self.width / 2 < i < jg_qr_st:
                will_x = i
                will_y = query_height + self.width / 2 - self.circle(self.width / 2, jg_qr_st - i)
                x.append(will_x)
                y.append(will_y)
            elif jg_qr_st <= i <= jg_qr_ed:
                will_x = i
                will_y = query_mar
                x.append(will_x)
                y.append(will_y)
            else:
                will_x = i
                will_y = query_height + self.width / 2 - self.circle(self.width / 2, i - jg_qr_ed)
                x.append(will_x)
                y.append(will_y)

        # p2->p1
        pos2_x = x[-1]
        pos2_y = y[-1]
        pos1_x = x[0]
        pos1_y = y[0]
        t = np.arange(0, 1.01, 0.005)
        dx = pos1_x - pos2_x
        dy = pos1_y - pos2_y
        p1, p2, p3, p4 = pos2_x, dx * ratio + pos2_x, -dx * ratio + pos1_x, pos1_x
        x1 = bezier3(p1, p2, p3, p4, t)
        p1, p2, p3, p4 = pos2_y, (1 - ratio) * dy + pos2_y, -(1 - ratio) * dy + pos1_y, pos1_y
        y1 = bezier3(p1, p2, p3, p4, t)
        x += x1
        y += y1

        return x, y

    def loose_compact_qry_ref(self, chr_plus_gap_length, length_df, length_df_chr_list, max_gm_gap_length):
        if self.gap_style == "loose":
            gap_length = (chr_plus_gap_length - length_df['length'].sum()) / (len(length_df_chr_list) + 1)
            start_list, end_list, chr_to_start = self.line_get_pos_list_compact(length_df, gap_length, chr_plus_gap_length, flag=False)
        else:
            gap_length = max_gm_gap_length
            start_list, end_list, chr_to_start = self.line_get_pos_list_compact(length_df, gap_length, chr_plus_gap_length, flag=True)
        # 20240818 only_for_chr_margin_collinearity_plot
        cap_judge = zip(length_df_chr_list, start_list, end_list)
        cap_judge_chr = {}
        for ch, start, end in cap_judge:
            start_x = start / chr_plus_gap_length + self.width / 2
            end_x = end / chr_plus_gap_length - self.width / 2
            cap_judge_chr[ch] = [start_x, end_x]
        return cap_judge_chr, chr_to_start, start_list, end_list

    def sub_run(self, collinearity, prefix, ref_length, query_length, total_length, query_height, ref_height, loop, last_loop, strip_chr_abbr, sp_chr_color_dict, max_genome_gap_length):
        ref_chr_list = ref_length['chr'].tolist()
        query_chr_list = query_length['chr'].tolist()
        chr_color_dict = {}
        cmap = self.get_my_cmap()

        i = 1
        for ch in ref_chr_list:
            color_dict = {'chr': cmap(round(i / len(ref_chr_list), 2))}
            chr_color_dict[ch] = color_dict
            i += 1
        # figure geometry info
        cap_judge_ref_chr, ref_chr_to_start, ref_start_list, ref_end_list = self.loose_compact_qry_ref(
            total_length, ref_length, ref_chr_list, max_genome_gap_length)

        # relative length 1
        label_x = -0.1
        label_y = ref_height + self.width / 2
        if self.italic:
            plt.text(label_x, label_y, prefix[0], ha="center", va="center", fontsize=self.species_name_font_size, color='black', fontstyle='italic')
        else:
            plt.text(label_x, label_y, prefix[0], ha="center", va="center", fontsize=self.species_name_font_size,
                     color='black')
        for i in range(len(ref_chr_list)):
            ref_start_x = ref_start_list[i] / total_length
            ref_end_x = ref_end_list[i] / total_length
            x, y = self.plot_line_chr(ref_start_x, ref_end_x, self.width, ref_height)

            plt.fill(x, y, facecolor = sp_chr_color_dict[prefix[0]], alpha=.7, edgecolor = sp_chr_color_dict[prefix[0]], zorder=2)
            label_x = (ref_start_x + ref_end_x) / 2
            label_y = ref_height + self.width / 2
            text_chr = ref_chr_list[i][len(prefix[0]):]
            for abbr in strip_chr_abbr:
                if text_chr.startswith(abbr):
                    text_chr = text_chr[len(abbr):]
            if not self.hide_chr:
                plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.chr_font_size, color='black', zorder=3)
            else:
                pass

        cap_judge_query_chr, query_chr_to_start, query_start_list, query_end_list = self.loose_compact_qry_ref(
            total_length, query_length, query_chr_list, max_genome_gap_length)

        if loop == last_loop:
            label_x = -0.1
            label_y = query_height + self.width / 2
            if self.italic:
                plt.text(label_x, label_y, prefix[1], ha="center", va="center", fontsize=self.species_name_font_size, color='black', fontstyle='italic')
            else:
                plt.text(label_x, label_y, prefix[1], ha="center", va="center", fontsize=self.species_name_font_size,
                         color='black')
            for i in range(len(query_chr_list)):
                query_start_x = query_start_list[i] / total_length
                query_end_x = query_end_list[i] / total_length
                x, y = self.plot_line_chr(query_start_x, query_end_x, self.width, query_height)

                plt.fill(x, y, facecolor = sp_chr_color_dict[prefix[1]], alpha=0.7, edgecolor=sp_chr_color_dict[prefix[1]], zorder=2)
                label_x = (query_start_x + query_end_x) / 2
                label_y = query_height + self.width / 2
                text_chr = query_chr_list[i][len(prefix[1]):]
                for abbr in strip_chr_abbr:
                    if text_chr.startswith(abbr):
                        text_chr = text_chr[len(abbr):]
                if not self.hide_chr:
                    plt.text(label_x, label_y, text_chr, ha="center", va="center", fontsize=self.chr_font_size, color='black', zorder=3)
                else:
                    pass

        chr_list = list(OrderedDict.fromkeys(query_chr_list + ref_chr_list))
        chr_to_start = {}
        chr_to_start.update(query_chr_to_start)
        chr_to_start.update(ref_chr_to_start)

        if self.actual_len:
            data, rf_blk_chr, qry_blk_chr, direction_list = self.read_proali_collinearity(prefix[1], prefix[0], collinearity, chr_list, chr_to_start)
        else:
            ratio_pair = self.ratio_pair[loop]
            data, rf_blk_chr, qry_blk_chr, direction_list = self.read_proali_collinearity_ratio(prefix[1], prefix[0], collinearity, chr_list, chr_to_start, ratio_pair)


        intra = []
        intra_chr_list = []
        intra_direction_list = []
        i = 0
        for block in data:
            if qry_blk_chr[i] == rf_blk_chr[i]:
                intra.append(block)
                intra_direction_list.append(direction_list[i])
                intra_chr_list.append(qry_blk_chr[i])
                i += 1
            else:
                color = chr_color_dict[rf_blk_chr[i]]['chr']
                pos1, pos2, pos3, pos4 = block[0], block[1], block[2], block[3]
                pos1_coord_x = pos1 / total_length
                pos2_coord_x = pos2 / total_length
                pos3_coord_x = pos3 / total_length
                pos4_coord_x = pos4 / total_length

                query_chr = qry_blk_chr[i]
                ref_chr = rf_blk_chr[i]
                judge_fake_query_start_x = cap_judge_query_chr[query_chr][0]
                judge_fake_query_end_x = cap_judge_query_chr[query_chr][1]
                judge_fake_ref_start_x = cap_judge_ref_chr[ref_chr][0]
                judge_fake_ref_end_x = cap_judge_ref_chr[ref_chr][1]
                x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height,
                                                     judge_fake_query_start_x, judge_fake_query_end_x, judge_fake_ref_start_x, judge_fake_ref_end_x)
                if self.color_style == "two_colors" and direction_list[i] == "+":
                    color = "#F0F0F0"
                    plt.fill(x, y, facecolor=color, alpha=0.7, zorder=1.5)
                elif self.color_style == "two_colors" and direction_list[i] == "-":
                    color = "#66AD56"
                    plt.fill(x, y, facecolor=color, alpha=0.7, zorder=1.8)
                else:
                    plt.fill(x, y, facecolor=color, alpha=0.7, zorder=1.8)
                i += 1
        i = 0
        for block in intra:
            color = chr_color_dict[intra_chr_list[i]]['chr']
            pos1, pos2, pos3, pos4 = block[0], block[1], block[2], block[3]
            pos1_coord_x = pos1 / total_length
            pos2_coord_x = pos2 / total_length
            pos3_coord_x = pos3 / total_length
            pos4_coord_x = pos4 / total_length

            query_chr = ref_chr = intra_chr_list[i]
            judge_fake_query_start_x = cap_judge_query_chr[query_chr][0]
            judge_fake_query_end_x = cap_judge_query_chr[query_chr][1]
            judge_fake_ref_start_x = cap_judge_ref_chr[ref_chr][0]
            judge_fake_ref_end_x = cap_judge_ref_chr[ref_chr][1]
            x, y = self.plot_collinearity_region(pos1_coord_x, pos2_coord_x, pos3_coord_x, pos4_coord_x, query_height, ref_height,
                                                 judge_fake_query_start_x, judge_fake_query_end_x, judge_fake_ref_start_x, judge_fake_ref_end_x)
            if self.color_style == "two_colors" and intra_direction_list[i] == "+":
                color = "#F0F0F0"
                plt.fill(x, y, facecolor=color, alpha=0.7, zorder=1.5)
            elif self.color_style == "two_colors" and intra_direction_list[i] == "-":
                color = "#66AD56"
                plt.fill(x, y, facecolor=color, alpha=0.7, zorder=1.8)
            else:
                plt.fill(x, y, facecolor=color, alpha=0.7, zorder=1.8)
            i += 1

    def run(self):
        sp_number = len(self.split_conf(self.collinearity, ",")) + 1
        self.determine_fig_par(sp_number)
        strip_chr_abbr = self.split_conf(self.remove_chromosome_prefix, ",")
        strip_collinearity_file_list = self.split_conf(self.collinearity, ",")

        # split prefix in order to pair
        prefix = self.prefix.split(",")
        strip_prefix = []
        new_prefix = []
        for prx in prefix:
            if len(prx) == 0:
                continue
            prx = prx.strip()
            strip_prefix.append(prx)
        for i in range(len(strip_collinearity_file_list)):
            new_prefix.append([strip_prefix[i], strip_prefix[i + 1]])

        sp_chr_color_dict = self.get_sp_chr_color_dict(strip_prefix)

        if self.actual_len:
            # split length get dataframe
            length_file_list = self.length_file.split(",")
            strip_length_file_list_df = []

            i = 0
            for le in length_file_list:
                le = le.strip()
                if len(le) == 0:
                    continue
                length = pd.read_csv(le, sep='\t', header=0)
                # print(length)
                length['chr'] = length['chr'].astype(str)
                length['chr'] = strip_prefix[i] + length['chr']
                length['length'] = length['length'].astype(int)
                strip_length_file_list_df.append(length)
                i += 1
        else:
            length_list = self.get_sp_length()
            ratio_list = self.get_expand_ratio(length_list)
            strip_length_file_list_df = self.split_length_adjust(strip_prefix, ratio_list)

            self.ratio_pair = []
            for i in range(len(strip_collinearity_file_list)):
                self.ratio_pair.append([ratio_list[i], ratio_list[i + 1]])

        new_length_file_list_df = []
        for i in range(len(strip_collinearity_file_list)):
            new_length_file_list_df.append([strip_length_file_list_df[i], strip_length_file_list_df[i+1]])

        # get maximum chr total length
        length_list = []
        for df in strip_length_file_list_df:
            length_list.append(df['length'].sum())
        total_chr_length = max(length_list)
        total_length = total_chr_length + total_chr_length / self.gap_ratio

        mx_chr_number = 10
        for df in strip_length_file_list_df:
            if df['length'].sum() == total_chr_length:
                mx_chr_number = len(df['chr'])
        max_gm_gap_length = total_length / self.gap_ratio / mx_chr_number

        zipped_three_pair = list(zip(strip_collinearity_file_list, new_prefix, new_length_file_list_df))
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'Bitstream Vera Serif',
                                      'Computer Modern Roman',
                                      'New Century Schoolbook', 'Century Schoolbook L', 'Utopia',
                                      'ITC Bookman', 'Bookman', 'Nimbus Roman No9 L', 'Times', 'Palatino', 'Charter',
                                      'serif']
        if self.figsize:
            fig_lt = self.figsize.split(",")
            fig, ax = plt.subplots(figsize=(float(fig_lt[0]), float(fig_lt[1])), facecolor='white')
        else:
            fig, ax = plt.subplots(figsize=(14, 14), facecolor='white')

        ax.set_aspect('equal')

        query_height = self.query_height
        ref_height = self.ref_height
        last_loop = len(zipped_three_pair) - 1
        i = 0
        for col, prefix, length in zipped_three_pair:
            self.sub_run(col, prefix, length[0], length[1], total_length, query_height, ref_height, i, last_loop, strip_chr_abbr, sp_chr_color_dict, max_gm_gap_length)
            query_height = query_height + self.height_gap

            ref_height = ref_height + self.height_gap
            i += 1

        plt.axis('off')
        # # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
        plt.savefig(self.savefig, dpi=int(self.dpi), bbox_inches='tight')
        print("produce", self.savefig, "success")
        sys.exit(0)
