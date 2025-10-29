#python select.block.py 47.Chr21_vs_CAAScds.Bgenome.anchors filtered.block.txt 47HCLC.final.4B.gff3 CSHCLC.gene.4B.gff3
import sys
import matplotlib.pyplot as plt
import numpy as np
import math


def bezier3(a, b, c, d, t):
    fvalue = []
    for i in t:
        y = a * (1-i)**3 + 3 * b * i * (1-i)**2 + 3 * c * (1-i) * i**2 + d * i**3
        fvalue.append(y)
    return fvalue

def plot_collinearity_region(pos1_x_, pos2_x_, pos3_x_, pos4_x_, pos1_y_, pos2_y_, pos3_y_, pos4_y_, _ratio):
    ratio = _ratio
    x_lt, y_lt = [], []

    # p1->p2(longest left -> right)
    x_lt.append(pos1_x_)
    y_lt.append(pos1_y_)

    x_lt.append(pos2_x_)
    y_lt.append(pos2_y_)

    # p2->p3
    x_lt.append(pos2_x_)
    y_lt.append(pos2_y_)
    t = np.arange(0, 1.01, 0.005)
    dx = pos3_x_ - pos2_x_
    dy = pos3_y_ - pos2_y_
    p1, p2, p3, p4 = pos2_x_, dx * ratio + pos2_x_, -dx * ratio + pos3_x_, pos3_x_
    x1 = bezier3(p1, p2, p3, p4, t)
    p1, p2, p3, p4 = pos2_y_, (1 - ratio) * dy + pos2_y_, -(1 - ratio) * dy + pos3_y_, pos3_y_
    y1 = bezier3(p1, p2, p3, p4, t)
    x_lt += x1
    y_lt += y1

    # p3->p4
    x_lt.append(pos3_x_)
    y_lt.append(pos3_y_)
    x_lt.append(pos4_x_)
    y_lt.append(pos4_y_)

    # p4->p1
    x_lt.append(pos4_x_)
    y_lt.append(pos4_y_)
    t = np.arange(0, 1.01, 0.005)
    dx = pos1_x_ - pos4_x_
    dy = pos1_y_ - pos4_y_
    p1, p2, p3, p4 = pos4_x_, dx * ratio + pos4_x_, -dx * ratio + pos1_x_, pos1_x_
    x1 = bezier3(p1, p2, p3, p4, t)
    p1, p2, p3, p4 = pos4_y_, (1 - ratio) * dy + pos4_y_, -(1 - ratio) * dy + pos1_y_, pos1_y_
    y1 = bezier3(p1, p2, p3, p4, t)
    x_lt += x1
    y_lt += y1

    return x_lt, y_lt

def get_rct_coord(x, y, length, width):
    x_inn = []
    y_inn = []
    x_inn.extend([x, x + length, x + length, x])
    y_inn.extend([y, y, y + width, y + width])
    return x_inn, y_inn

def plot_multiple_rct(coord_list):
    for rct_x, rct_y, rct_length, rct_width, color, alpha_sub, zorder_sub in coord_list:
        x_sub, y_sub = get_rct_coord(rct_x, rct_y, rct_length, rct_width)
        plt.fill(x_sub, y_sub, facecolor=color, alpha=alpha_sub, zorder=zorder_sub)

def plot_synteny(x_list, y_list, facecolor_sub="#F0F0F0", alpha_sub=0.51, zorder_sub=2, ratio=0.5):
    pos1_x, pos2_x, pos3_x, pos4_x = x_list
    pos1_y, pos2_y, pos3_y, pos4_y = y_list
    x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y, ratio)
    plt.fill(x_, y_, facecolor=facecolor_sub, alpha=alpha_sub, zorder=zorder_sub)

plt.rcParams['font.family'] = "Times New Roman"
fig, ax = plt.subplots(figsize=(14, 8), facecolor='white')
ax.set_aspect('equal')

anchors = sys.argv[1]
filtered_block = sys.argv[2]
output = open(filtered_block, "w")

query_start_mx = 250 * 1000000
query_end_mx = 400 * 1000000

count = 0
loop_ratio = 0
with open(anchors) as f:
    next(f)
    next(f)
    flag = False
    for line in f:
        if line.startswith("#block begin"):
            # enter block
            flag = True
            block = []
            ref_coor = []
            output_block = []
            continue
        if line.startswith("#block end"):
            # off block
            block_start = min(block)
            block_end = max(block)
            ref_start = min(ref_coor)
            ref_end = max(ref_coor)
            if query_start_mx < block_start < block_end < query_end_mx:
                # output.write("#begin" + "\n")

                outer_x_list = [int(ref_coor[0]) * (10**-6), int(ref_coor[-1]) * (10**-6), int(block[-1]) * (10**-6), int(block[0]) * (10**-6)]
                outer_y_list = [0, 0, 30, 30]
                if ref_chr == query_chr == "Chr4B":
                    # print("query_start: ", int(block[0]) * (10 ** -6))
                    # print("query_end: ", int(block[-1]) * (10 ** -6))
                    # print("ref_start: ", int(ref_coor[0]) * (10 ** -6))
                    # print("ref_end: ", int(ref_coor[-1]) * (10 ** -6))
                    if block_direction == "+":
                        plot_synteny(outer_x_list, outer_y_list, facecolor_sub="#F0F0F0", alpha_sub=1, zorder_sub=2, ratio=0.316)
                        pass
                    else:
                        if count == 0:
                            loop_ratio = 0.318
                        elif count ==1:
                            loop_ratio = 0.318
                        else:
                            loop_ratio = 0.318
                        count += 1
                        plot_synteny(outer_x_list, outer_y_list, facecolor_sub="#66AD56", alpha_sub=0.51, zorder_sub=3, ratio=loop_ratio)
                print()
                # output.writelines(output_block)
                # output.write("#end" + "\n")
            flag = False
            continue
        if flag:
            output_block.append(line)
            row_list = line.split()
            query_start = row_list[4]
            query_end = row_list[5]
            ref_chr = row_list[0]
            query_chr = row_list[3]
            block_direction = row_list[6]
            block.append(int(query_start))
            block.append(int(query_end))
            ref_coor.append(int(row_list[1])+15042850)
            ref_coor.append(int(row_list[2])+15042850)


outer_rct_width = 1
outer_coord_list = [(251.848292, 30, 137.035727, outer_rct_width, "#9C80EE", 1, 2),
                    (251.848292, -1,133.167494 , outer_rct_width, "#F8B62D", 1, 2)]
plot_multiple_rct(outer_coord_list)


ref_gff = sys.argv[3]
query_gff = sys.argv[4]

with open(ref_gff) as f:
    for ln in f:
        row_lt = ln.split()
        gn_st = int(row_lt[3])
        gn_end = int(row_lt[4])
        sec_column = row_lt[2]
        eight_column = row_lt[8].split("ID=")[1].strip()
        print(gn_st)
        print(gn_end)
        print(sec_column)
        if 251848292 - 15042850 < gn_st and gn_end < 251848292 + 133167494 - 15042850:
            # print("aaaaaaaaaaaaaaaa")
            if sec_column == "gene":
                # print("bbbbbbbbbbbbbbbbb")
                if not eight_column.endswith("LC"):
                    # print("assssssssssssss")
                    outer_rct_width = 1
                    outer_coord_list = [((gn_st+15042850)*(10**-6), -1, (gn_end-gn_st)*(10**-6), outer_rct_width, "black", 1, 2)]
                    plot_multiple_rct(outer_coord_list)

with open(query_gff) as f:
    for ln in f:
        row_lt = ln.split()
        gn_st = int(row_lt[3])
        gn_end = int(row_lt[4])
        sec_column = row_lt[2]
        eight_column = row_lt[8].split("ID=")[1].split(";")[0]
        print(gn_st)
        print(gn_end)
        print(sec_column)
        if 251848292 < gn_st and gn_end < 137035727 + 251848292:
            # print("aaaaaaaaaaaaaaaa")
            if sec_column == "gene":
                # print("bbbbbbbbbbbbbbbbb")
                if eight_column.endswith("HC"):
                    # print("assssssssssssss")
                    outer_rct_width = 1
                    outer_coord_list = [((gn_st)*(10**-6), 30, (gn_end-gn_st)*(10**-6), outer_rct_width, "black", 1, 2)]
                    plot_multiple_rct(outer_coord_list)
plt.axis('off')
# # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
plt.savefig("./synteny.tandem.pdf", dpi=300)
