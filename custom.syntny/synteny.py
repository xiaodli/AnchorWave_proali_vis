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

def get_rct_coord(x, y):
    x_inn = []
    y_inn = []
    x_inn.extend([x, x + 591, x + 591, x])
    y_inn.extend([y, y, y + 5, y + 5])
    x = 0
    y = y + 5+150
    return x_inn, y_inn, x, y

def get_longest_rct_coord(x, y, length=300):
    x_inn = []
    y_inn = []
    x_inn.extend([x, x + length, x + length, x])
    y_inn.extend([y, y, y + 5, y + 5])
    x = 0
    y = y + 5+150
    return x_inn, y_inn, x, y

def plot_longest_gene():
    x_, y_, _, _ = get_longest_rct_coord(700, 930, 290)
    plt.fill(x_, y_, facecolor="blue", alpha=.7, zorder=2)

    x_, y_, _, _ = get_longest_rct_coord(990, 930, 2)
    plt.fill(x_, y_, facecolor="red", alpha=1, zorder=2)

    x_, y_, _, _ = get_longest_rct_coord(992, 930, 161)
    plt.fill(x_, y_, facecolor="green", alpha=.7, zorder=2)

    x_, y_, _, _ = get_longest_rct_coord(1153, 930, 50)
    plt.fill(x_, y_, facecolor="purple", alpha=.7, zorder=2)

    x_, y_, _, _ = get_longest_rct_coord(1203, 930, 30)
    plt.fill(x_, y_, facecolor="white", alpha=.7, zorder=2)

    x_, y_, _, _ = get_longest_rct_coord(1233, 930, 50)
    plt.fill(x_, y_, facecolor="purple", alpha=.7, zorder=2)

    x_, y_, _, _ = get_longest_rct_coord(1283, 930, 138)
    plt.fill(x_, y_, facecolor="pink", alpha=.7, zorder=2)

def plot_top_gene():
    # plot top gene
    x_, y_, _, _ = get_rct_coord(300, 1185)
    plt.fill(x_, y_, facecolor="red", alpha=.7, zorder=2)

def plot_three_gene():
    three_bottom_x = 765
    three_bottom_y = 465
    for i in range(3):
        x_, y_, _, _ = get_longest_rct_coord(three_bottom_x, three_bottom_y, 290)
        plt.fill(x_, y_, facecolor="blue", alpha=.7, zorder=2)

        x_, y_, _, _ = get_longest_rct_coord(three_bottom_x+290, three_bottom_y, 2)
        plt.fill(x_, y_, facecolor="red", alpha=1, zorder=2)

        x_, y_, _, _ = get_longest_rct_coord(three_bottom_x+292, three_bottom_y, 161)
        plt.fill(x_, y_, facecolor="green", alpha=.7, zorder=2)

        x_, y_, _, _ = get_longest_rct_coord(three_bottom_x+453, three_bottom_y, 138)
        plt.fill(x_, y_, facecolor="pink", alpha=.7, zorder=2)
        three_bottom_y += 155

def plot_seven_gene():
    left_bottom_x = 0
    left_bottom_y = 0
    for i in range(7):
        x_, y_, left_bottom_x, left_bottom_y = get_rct_coord(left_bottom_x, left_bottom_y)
        plt.fill(x_, y_, facecolor="red", alpha=.7, zorder=2)

def plot_seven_synteny():
    gap_length = 150
    bottom_y = 5
    for i in range(6):
        pos1_x, pos2_x, pos3_x, pos4_x = 0, 591, 591, 0
        pos1_y, pos2_y, pos3_y, pos4_y = bottom_y, bottom_y, bottom_y+gap_length, bottom_y+gap_length
        x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
        plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.51, zorder=2)
        bottom_y = bottom_y + 5 + gap_length
    pos1_x, pos2_x, pos3_x, pos4_x = 0, 591, 891, 300
    pos1_y, pos2_y, pos3_y, pos4_y = 935, 935, 1185, 1185
    x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
    plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.51, zorder=1)

def plot_three_synteny():
    gap_length = 150
    bottom_y = 470
    for i in range(2):
        pos1_x, pos2_x, pos3_x, pos4_x = 765, 1356, 1356, 765
        pos1_y, pos2_y, pos3_y, pos4_y = bottom_y, bottom_y, bottom_y+gap_length, bottom_y+gap_length
        x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
        plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.51, zorder=2)
        bottom_y = bottom_y + 5 + gap_length

    bottom_y = 780
    pos1_x, pos2_x, pos3_x, pos4_x = 765, 1218, 1153, 700
    pos1_y, pos2_y, pos3_y, pos4_y = bottom_y, bottom_y, bottom_y+gap_length, bottom_y+gap_length
    x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
    plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.81, zorder=2)

    pos1_x, pos2_x, pos3_x, pos4_x = 1218, 1356, 1421, 1283
    pos1_y, pos2_y, pos3_y, pos4_y = bottom_y, bottom_y, bottom_y + gap_length, bottom_y + gap_length
    x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
    plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.81, zorder=2)

######################################################################
plt.rcParams['font.family'] = "Times New Roman"
fig, ax = plt.subplots(figsize=(14, 8), facecolor='white')
ax.set_aspect('equal')

## plot gene
plot_longest_gene()
plot_three_gene()
plot_seven_gene()
plot_top_gene()

## plot synteny left
pos1_x, pos2_x, pos3_x, pos4_x = 700, 990, 590, 300
pos1_y, pos2_y, pos3_y, pos4_y = 935, 935, 1185, 1185
x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.81, zorder=2)

pos1_x, pos2_x, pos3_x, pos4_x = 992, 1153, 751, 590
pos1_y, pos2_y, pos3_y, pos4_y = 935, 935, 1185, 1185
x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.81, zorder=2)

pos1_x, pos2_x, pos3_x, pos4_x = 1283, 1421, 891, 751,
pos1_y, pos2_y, pos3_y, pos4_y = 935, 935, 1185, 1185
x_, y_ = plot_collinearity_region(pos1_x, pos2_x, pos3_x, pos4_x, pos1_y, pos2_y, pos3_y, pos4_y)
plt.fill(x_, y_, facecolor="#F0F0F0", alpha=0.81, zorder=2)

## plot normal synteny
plot_seven_synteny()

## plot three synteny
plot_three_synteny()

y1 = 932.5 - math.sqrt(3)*(1198-1203)
y2 = 932.5 - math.sqrt(3)*(1208-1203)
plt.plot([1198,1208], [y1,y2], color="b")

y1 = 932.5 - math.sqrt(3)*(1228-1233)
y2 = 932.5 - math.sqrt(3)*(1238-1233)
plt.plot([1228,1238], [y1,y2], color="b")

plt.axis('off')
# # plt.subplots_adjust(0.1, 0.1, 0.9, 0.9)
plt.savefig("/home/dell/Desktop/other_people/synteny_wheat_line/test.png", dpi=300)
