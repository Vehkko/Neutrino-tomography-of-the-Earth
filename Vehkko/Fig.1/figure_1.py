"""
cosθ 的范围取 [-1.0, 0.1], 共 60 个 bin
能量范围: [400, 20000] Gev
数据文件: data/observed_events
error bar 的处理: 区间内的事件数开根号作为绝对误差

"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# 字体设置
plt.rcParams["font.family"] = "Times New Roman"

# 数据文件
filename = "./data/observed_events.dat"

# 能量阈值
E_rec_inf, E_rec_sup = 400, 20000
E_rec0, E_rec1 = 1500, 2500


def read_energy_and_cosine(file_path):
    """从数据文件读取能量和天顶角余弦值"""
    energies, cos_zenith = [], []

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue  # 跳过注释行

            parts = line.split()
            if len(parts) == 2:
                try:
                    energies.append(float(parts[0]))  # 读取能量
                    cos_zenith.append(np.cos(float(parts[1])))  # 读取天顶角余弦
                except ValueError:
                    print(f"跳过无法解析的行: {line.strip()}")

    return np.array(energies), np.array(cos_zenith)


# 读取数据
energies, cos_zenith = read_energy_and_cosine(filename)

# 筛选符合能量范围的数据
mask = (E_rec_inf <= energies) & (energies <= E_rec_sup)
energies, cos_zenith = energies[mask], cos_zenith[mask]

print("Number of events:", len(cos_zenith))

# 定义 bins
bins = np.linspace(-1, 0.1, 60)


# 计算不同能量区间的分布
def compute_histogram(energies, cos_zenith, bins, E_min, E_max):
    mask = (E_min <= energies) & (energies <= E_max)
    return np.histogram(cos_zenith[mask], bins=bins)


hist1_raw, bin_edges1 = compute_histogram(
    energies, cos_zenith, bins, E_rec_inf, E_rec_sup
)
hist2_raw, bin_edges2 = compute_histogram(energies, cos_zenith, bins, E_rec0, E_rec_sup)
hist3_raw, bin_edges3 = compute_histogram(energies, cos_zenith, bins, E_rec1, E_rec_sup)

# 绘制图像
# plt.figure(figsize=(8, 6))
fig, ax = plt.subplots(figsize=(4.8, 6.3))

# 背景色区域
bg_regions = [
    (-1, -0.98, "#ffb2b2"),
    (-0.98, -0.82, "#d9b9ab"),
    (-0.82, 0.0, "#e5cccc"),
    (0.0, 0.2, "#e3e3ff"),
]

for x_min, x_max, color in bg_regions:
    plt.axvspan(x_min, x_max, facecolor=color, alpha=0.75)

# 壳层分界线
plt.axvline(-0.98, color="#d0715f", linestyle="-", linewidth=1)
plt.axvline(-0.82, color="#b67368", linestyle="-", linewidth=1)
plt.axvline(0.0, color="#b58ead", linestyle="-", linewidth=1)

# 壳层标签
plt.text(
    -0.976,
    3.00,
    "InnerCore",
    color="#f81a17",
    fontsize=12,
    ha="center",
    va="center",
    rotation=90,
)
plt.text(
    -0.84,
    3.00,
    "OuterCore",
    color="#ae4643",
    fontsize=12,
    ha="center",
    va="center",
    rotation=90,
)
plt.text(
    -0.02,
    3.05,
    "Mantle",
    color="#9a3636",
    fontsize=12,
    ha="center",
    va="center",
    rotation=90,
)
plt.text(
    0.08,
    2.96,
    "Atmosphere",
    color="#5b5bb9",
    fontsize=12,
    ha="center",
    va="center",
    rotation=90,
)


def plot_with_error_bars(
    bin_edges, hist_raw, color, lable_pos_x, lable_pos_y, lable_text
):
    x = bin_edges[:-1] + np.diff(bin_edges) / 2  # 计算 bin 的中心
    hist_log = np.log10(hist_raw, where=hist_raw > 0)  # 计算 log10 值

    # 计算泊松误差棒
    sigma_plus = np.log10(1 + 1 / np.sqrt(hist_raw))
    sigma_minus = -np.log10(1 - 1 / np.sqrt(hist_raw))

    # 矩形宽度
    width = np.diff(bin_edges)[0]

    # 用于绘制竖直线
    y_last = 0

    for xi, yi, s_m, s_p in zip(x, hist_log, sigma_minus, sigma_plus):
        # 绘制误差棒的小方框
        rect = Rectangle(
            (xi - 0.5 * width, yi - s_m), width, s_p + s_m, color="#7f65cc"
        )
        ax.add_patch(rect)

        # 绘制小方框之间的竖直线
        plt.plot(
            [xi - 0.5 * width, xi - 0.5 * width],
            ([y_last, yi] if y_last < yi else [yi, y_last]),
            color="#7171ff",
            linewidth=1.2,
        )
        y_last = yi

        # 绘制误差棒左右两端的小短线
        plt.plot(
            [xi - 0.5 * width, xi - 0.5 * width],
            [yi - s_m, yi + s_p],
            color="#7171ff",
            linewidth=0.8,
        )
        plt.plot(
            [xi + 0.5 * width, xi + 0.5 * width],
            [yi - s_m, yi + s_p],
            color="#7171ff",
            linewidth=1.2,
        )

        # 绘制误差棒上下两端的小短线
        plt.plot(
            [xi - 0.5 * width - 0.0001, xi + 0.5 * width + 0.0001],
            [yi - s_m, yi - s_m],
            color,
            linewidth=0.8,
        )
        plt.plot(
            [xi - 0.5 * width - 0.0001, xi + 0.5 * width + 0.0001],
            [yi + s_p, yi + s_p],
            color,
            linewidth=0.8,
        )

    # 绘制散点和误差棒
    plt.errorbar(
        x,
        hist_log,
        yerr=[sigma_minus, sigma_plus],
        fmt="o",
        color=color,
        markersize=2,
        elinewidth=1,
        capsize=0,
    )

    # 在曲线旁边添加标签
    plt.text(
        lable_pos_x,
        lable_pos_y,
        lable_text,
        color="black",
        fontsize=10,
        fontweight="bold",
        ha="center",
        va="center",
        rotation=10,
    )


# 绘制三组散点图
plot_with_error_bars(bin_edges1, hist1_raw, "blue", -0.54, 2.55, r"$Full~sample$")
plot_with_error_bars(
    bin_edges2, hist2_raw, "blue", -0.54, 2.10, r"$E_{\mu}^{rec} > 1.5 TeV$"
)
plot_with_error_bars(
    bin_edges3, hist3_raw, "blue", -0.54, 1.62, r"$E_{\mu}^{rec} > 2.5 TeV$"
)

# 图例
# plt.legend(loc="upper right")

# 坐标轴设置
plt.xlabel(r"$\cos \theta_z^{rec}$", fontsize=12)
plt.ylabel(r"$\lg(N_{events})$", fontsize=12)
plt.ylim(0.77, 3.19)
plt.xlim(-1.0, 0.1)

# 保存图像
plt.savefig("figure/FIG.1.b.pdf")
plt.savefig("figure/FIG.1.b.png", dpi=300)

# 显示图像
plt.show()
