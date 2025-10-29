import matplotlib.pyplot as plt
import numpy as np
import math


def bezier3(a, b, c, d, t):
    fvalue = []
    for i in t:
        y = a * (1-i)**3 + 3 * b * i * (1-i)**2 + 3 * c * (1-i) * i**2 + d * i**3
        fvalue.append(y)
    return fvalue

def plot_collinearity_region(pos1_x_, pos2_x_, pos3_x_, pos4_x_, pos1_y_, pos2_y_, pos3_y_, pos4_y_):
    ratio = 0.316
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

def plot_synteny(x_list, y_list, facecolor_sub="#F0F0F0", alpha_sub=0.51, zorder_sub=2):
    pos1_x, pos2_x, pos3_x, pos4_x = x_list
    pos1_y, pos2_y, pos3_y, pos4_y = y_list
    x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
    plt.fill(x_, y_, facecolor=facecolor_sub, alpha=alpha_sub, zorder=zorder_sub)

if __name__ == "__main__":
    ######################################################################
    plt.rcParams['font.family'] = "Times New Roman"
    fig, ax = plt.subplots(figsize=(14, 8), facecolor='white')
    ax.set_aspect('equal')

    ## plot top left gene and bottom left gene
    outer_rct_width = 100
    outer_coord_list = [(51559, 2100, 156, outer_rct_width, "#9C80EE", 1, 2),
                        (51715, 2100, 2512, outer_rct_width, "#F8B62D", 1, 2),
                        (54227, 2100, 309, outer_rct_width, "#9C80EE", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_coord_list = [(51559, 0, 156, outer_rct_width, "#9C80EE", 1, 2),
                        (51715, 0, 2512, outer_rct_width, "#F8B62D", 1, 2),
                        (54227, 0, 309, outer_rct_width, "#9C80EE", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_x_list = [51559, 54536, 54536, 51559]
    outer_y_list = [100, 100, 2100, 2100]
    plot_synteny(outer_x_list, outer_y_list, facecolor_sub="#F0F0F0", alpha_sub=1, zorder_sub=2)

    ## plot inter-gene region
    outer_rct_width = 50
    outer_coord_list = [(54536, 2125, 1000, outer_rct_width, "#C9CACA", 1, 2),
                        (55736, 2125, 1000, outer_rct_width, "#C9CACA", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_coord_list = [(54536, 25, 1000, outer_rct_width, "#C9CACA", 1, 2),
                        (55736, 25, 1000, outer_rct_width, "#C9CACA", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_x_list = [54536, 56736, 56736, 54536]
    outer_y_list = [100, 100, 2100, 2100]
    plot_synteny(outer_x_list, outer_y_list, facecolor_sub="#F0F0F0", alpha_sub=1, zorder_sub=2)

    y1 = 2150 - math.sqrt(3) * (55486 - 55536)
    y2 = 2150 - math.sqrt(3) * (55586 - 55536)
    plt.plot([55486, 55586], [y1, y2], color="b")

    y1 = 2150 - math.sqrt(3) * (55686 - 55736)
    y2 = 2150 - math.sqrt(3) * (55786 - 55736)
    plt.plot([55686, 55786], [y1, y2], color="b")

    y1 = 50 - math.sqrt(3) * (55486 - 55536)
    y2 = 50 - math.sqrt(3) * (55586 - 55536)
    plt.plot([55486, 55586], [y1, y2], color="b")

    y1 = 50 - math.sqrt(3) * (55686 - 55736)
    y2 = 50 - math.sqrt(3) * (55786 - 55736)
    plt.plot([55686, 55786], [y1, y2], color="b")

    ## plot top right gene and bottom right gene (UTR CDS)
    outer_rct_width = 100
    offset_inter = 3836
    outer_coord_list = [(60572-offset_inter, 2100, 154, outer_rct_width, "#9C80EE", 1, 2),
                        (60726-offset_inter, 2100, 2449, outer_rct_width, "#F8B62D", 1, 2),
                        (63175-offset_inter, 2100, 240, outer_rct_width, "#9C80EE", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_coord_list = [(60572-offset_inter, 0, 154, outer_rct_width, "#9C80EE", 1, 2),
                        (60726-offset_inter, 0, 2242, outer_rct_width, "#F8B62D", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_x_list = [60572-offset_inter, 62968-offset_inter, 62968-offset_inter, 60572-offset_inter]
    outer_y_list = [100, 100, 2100, 2100]
    plot_synteny(outer_x_list, outer_y_list, facecolor_sub="#F0F0F0", alpha_sub=1, zorder_sub=2)

    ## plot inter-gene region (right)
    outer_rct_width = 50
    outer_coord_list = [(63415-offset_inter, 2125, 1653, outer_rct_width, "#C9CACA", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_coord_list = [(62968-offset_inter, 25, 500, outer_rct_width, "#C9CACA", 1, 2)]
    plot_multiple_rct(outer_coord_list)

    outer_x_list = [62968-offset_inter, 63468-offset_inter, 65068-offset_inter, 64568-offset_inter]
    outer_y_list = [100, 100, 2100, 2100]
    plot_synteny(outer_x_list, outer_y_list, facecolor_sub="#F0F0F0", alpha_sub=1, zorder_sub=2)

    plt.axis('off')
    # # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
    plt.savefig("/home/dell/Desktop/other_people/synteny_wheat_line/synteny.tandem.svg", dpi=300)
