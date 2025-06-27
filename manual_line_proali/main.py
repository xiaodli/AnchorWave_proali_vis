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
import wheat_line_proali_pangenome
from pathlib import Path
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
    config_par.read(parameter.conf)
    file_empty(parameter.conf)
    wheat_line_proali_pangenome.Line(config_par).run()

parser = argparse.ArgumentParser(description='Conduct strand and WGD aware syntenic gene identification for a pair of genomes using the longest path algorithm implemented in AnchorWave.', prog="quota_Anchor")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.0.1')
subparsers = parser.add_subparsers(title='Gene collinearity analysis', dest='analysis')
# collinearity AnchorWave proali anchors plot
parser_sub_line_proali = subparsers.add_parser('line_proali', help='Anchors file from AnchorWave proali to visualization(line style).')
parser_sub_line_proali.set_defaults(func=run_line)
parser_sub_line_proali.add_argument('-c', '--conf', dest='conf', help="Configure file.", metavar="")
args = parser.parse_args()
args.func(args)