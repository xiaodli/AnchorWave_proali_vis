# -*- encoding: utf-8 -*-
'''
@File    :   main.py
@Time    :   2025/03/23 21:22:45
@Author  :   xiaodong li
@Version :   1.0.0
@Contact :   xiaodongli2405@gmail.com
'''


import configparser
import argparse
import line_proali_pangenome
import os
import sys


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

def run_line(parameter):
    config_par = configparser.ConfigParser()
    if parameter.conf is not None:
        file_empty(parameter.conf)
        config_par.read(parameter.conf)
    line_proali_pangenome.Line(config_par, parameter).run()

parser = argparse.ArgumentParser(description='AnchorWave proali subcommand anchors plot.', prog="quota_Anchor")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')
subparsers = parser.add_subparsers(title='Proali anchors visualization', dest='analysis')
# AnchorWave proali anchors plot
parser_sub_line_proali = subparsers.add_parser('line_proali', help='Anchors file from AnchorWave proali to visualization(line style).')
parser_sub_line_proali.set_defaults(func=run_line)
parser_sub_line_proali.add_argument('-c', '--conf', dest='conf', help="Configure file.", metavar="")
parser_sub_line_proali.add_argument('-i', '--collinearity', dest='collinearity', help="Collinearity file(e.g. file1,file2,file3)(Separator: ',').", metavar="")
parser_sub_line_proali.add_argument('-o', '--savefig', dest='savefig', help="Specify a file name to save figure.", metavar="")
parser_sub_line_proali.add_argument('-l', '--length_file', dest='length_file', help="Species length file list(e.g. file1,file2,file3)(Separator: ',').", metavar="")
parser_sub_line_proali.add_argument('-n', '--prefix', dest='prefix', help="Species name list(e.g. name1,name2,name3)(Separator: ',').", metavar="")
parser_sub_line_proali.add_argument('-rm', '--remove_chromosome_prefix', dest='remove_chromosome_prefix', help="Remove chromosome prefix to plot(e.g. chr,Chr,CHR)(Separator: ',').", metavar="")
parser_sub_line_proali.add_argument('-cf', '--chr_font_size', dest='chr_font_size', help="Chromosome name font size(defaults: 7).", type=int, metavar="")
parser_sub_line_proali.add_argument('-sf', '--species_name_font_size', dest='species_name_font_size', help="Species name font size(defaults: 7).", type=int, metavar="")
parser_sub_line_proali.add_argument('-hc', '--hide_chr', dest='hide_chr', help="Hide chromosome name in the figure.", action='store_true')
parser_sub_line_proali.add_argument('-it', '--italic', dest='italic', help="Species name italic in the figure.", action='store_true')
parser_sub_line_proali.add_argument('-sc', '--species_color', dest='sp_chr_color_comma_sep', help="Optional, species chromosome color(comma separated)", metavar="")
parser_sub_line_proali.add_argument('-fs', '--figsize', dest='figsize', help="Figure size(defaults: 14,14).", metavar="")
parser_sub_line_proali.add_argument('-cs', '--color_style', dest='color_style', help="Optional, block color style(rainbow, husl, four_colors, two_colors).", type=str, choices=["rainbow", "husl", "four_colors", "two_colors"], metavar="")
parser_sub_line_proali.add_argument('-al', '--actual_len', dest='actual_len', help="Optional, use actual chromosome length in the figure.", action='store_true')
parser_sub_line_proali.add_argument('-gs', '--gap_style', dest='gap_style', help="Typing compact or loose(default: loose).", type=str, choices=["loose", "compact"],  metavar="")
args = parser.parse_args()
args.func(args)
