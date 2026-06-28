import os
import re
import numpy as np
import matplotlib.pyplot as plt

REACTION = "ArBe_qgg"
FILENAME = "ArBe_qgg_v3.dat"


def read_data_file(filename):
    with open(filename, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    print(f"Прочитан файл: {filename}")

    n_groups = int(lines[0])
    counts = list(map(int, lines[1].split()))

    if len(counts) != n_groups:
        raise ValueError("Количество групп не совпадает с числом элементов во 2-й строке")

    groups = []
    current_line = 2

    for g in range(n_groups):
        group_name = f"Group_{g+1}"

        if current_line < len(lines) and lines[current_line].startswith("#"):
            header = lines[current_line]
            current_line += 1

            m = re.search(r'#\d+\s+([A-Za-z]+)', header)
            if m:
                group_name = m.group(1)

        points = []

        for _ in range(counts[g]):
            if current_line >= len(lines):
                raise ValueError("Файл закончился раньше, чем ожидалось")

            parts = lines[current_line].split()
            current_line += 1

            if len(parts) < 6:
                raise ValueError(f"Некорректная строка данных: {lines[current_line-1]}")

            label = parts[1]          # 2-й столбец
            x = float(parts[2])       # 3-й столбец: Qgg
            n = float(parts[3])       # 4-й столбец: N
            y = float(parts[4])       # 5-й столбец: ln(N)
            yerr = float(parts[5])    # 6-й столбец: 1/sqrt(N)

            points.append((x, y, yerr, n, label))

        groups.append({
            "name": group_name,
            "points": points
        })

    return groups

def read_params_file(filename, n_groups):
    a = None
    b_dict = {}

    print(f"Прочитан файл: {filename}")

    with open(filename, "r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            m_a = re.fullmatch(r"a\s*=\s*(-?[0-9]*\.?[0-9]+)", line, re.IGNORECASE)
            if m_a:
                a = float(m_a.group(1))
                continue

            m_b = re.fullmatch(r"b\((\d+)\)\s*=\s*(-?[0-9]*\.?[0-9]+)", line, re.IGNORECASE)
            if m_b:
                idx = int(m_b.group(1))
                val = float(m_b.group(2))
                b_dict[idx] = val
                continue

    if a is None:
        raise ValueError("В файле параметров не найден коэффициент a")

    b = []
    for j in range(1, n_groups + 1):
        if j not in b_dict:
            raise ValueError(f"В файле параметров не найден коэффициент b({j})")
        b.append(b_dict[j])

    return a, b


def plot_groups(groups, a, b, output_dir="plots"):
    os.makedirs(output_dir, exist_ok=True)
    colors = plt.cm.tab10.colors

    for j, group in enumerate(groups):
        name = group["name"]
        points = group["points"]
        bj = b[j]
        color = colors[j % len(colors)]

        x_data = np.array([p[0] for p in points])
        y_data = np.array([p[1] for p in points])
        y_err = np.array([p[2] for p in points])
        labels = [p[4] for p in points]

        x_min = min(x_data.min(), 0) - 2
        x_max = max(x_data.max(), 0) + 2
        x_line = np.linspace(x_min, x_max, 200)
        y_line = a * x_line + bj

        plt.figure(figsize=(8, 5))

        plt.errorbar(
            x_data,
            y_data,
            yerr=y_err,
            fmt='s',
            color=color,
            ecolor=color,
            elinewidth=1,
            capsize=3,
            markersize=5,
            label='Экспериментальные данные'
        )

        for x, y, label in zip(x_data, y_data, labels):
            plt.annotate(
                str(label),
                (x, y),
                textcoords="offset points",
                xytext=(5, 5),
                ha='left',
                fontsize=9,
                color='black'
            )

        plt.plot(
            x_line,
            y_line,
            color=color,
            linewidth=1.5,
            label=f'ln(N) = {a:.3f} * Qgg + {bj:.3f}'
        )

        plt.title(name, fontsize=18, fontweight='bold')
        plt.xlabel('Qgg (МэВ)', fontsize=14, fontweight='bold', loc='right')
        plt.ylabel('ln(N)', fontsize=14, fontweight='bold', rotation=0, loc='top')
        plt.legend(frameon=False, loc='upper left')
        plt.tight_layout()

        out_path = os.path.join(output_dir, f"{j+1:02d}_{name}.png")
        plt.savefig(out_path, dpi=150)
        plt.close()

        print(f"Сохранен график: {out_path}")

def plot_all_groups(groups, a, b, output_dir="plots"):
    os.makedirs(output_dir, exist_ok=True)
    colors = plt.cm.tab10.colors

    plt.figure(figsize=(10, 7))

    global_xmin = float("inf")
    global_xmax = float("-inf")

    for j, group in enumerate(groups):
        name        = group["name"]
        points      = group["points"]
        bj          = b[j]
        color       = colors[j % len(colors)]

        x_data = np.array([p[0] for p in points])
        y_data = np.array([p[1] for p in points])
        y_err  = np.array([p[2] for p in points])
        labels = [p[4] for p in points]

        global_xmin = min(global_xmin, x_data.min())
        global_xmax = max(global_xmax, x_data.max())

        plt.errorbar(
            x_data,
            y_data,
            yerr=y_err,
            fmt='s',
            color=color,
            ecolor=color,
            elinewidth=1,
            capsize=3,
            markersize=5,
            linestyle='none',
            label=f'{name}'
        )

        for x, y, label in zip(x_data, y_data, labels):
            plt.annotate(
                str(label),
                (x, y),
                textcoords="offset points",
                xytext=(4, 4),
                ha='left',
                fontsize=10,
                color='black'
            )

    x_line = np.linspace(min(global_xmin, 0) - 2, max(global_xmax, 0) + 2, 300)

    for j, group in enumerate(groups):
        name = group["name"]
        bj = b[j]
        color = colors[j % len(colors)]
        y_line = a * x_line + bj

        plt.plot(
            x_line,
            y_line,
            color=color,
            linewidth=1.5,
            label=f'{name}: y={a:.3f}x+{bj:.3f}'
        )

    plt.xlabel('Qgg (МэВ)', fontsize=14, fontweight='bold')
    plt.ylabel('ln(N)', fontsize=14, fontweight='bold')
    plt.legend(frameon=False, fontsize=9, ncol=2)
    plt.tight_layout()

    out_path = os.path.join(output_dir, REACTION + "_all_groups.png")
    plt.savefig(out_path, dpi=150)
    plt.close()

    print(f"Сохранен общий график: {out_path}")

def main():
    data_file = "../input/" + FILENAME
    params_file = "../output/" + REACTION + "_out.dat"
    outdir = REACTION + "_plots"

    groups = read_data_file(data_file)
    a, b = read_params_file(params_file, len(groups))

    print(f"a = {a}")
    for i, bj in enumerate(b, start=1):
        print(f"b{i} = {bj}")

    plot_groups(groups, a, b, output_dir=outdir)
    plot_all_groups(groups, a, b, output_dir=outdir)


if __name__ == "__main__":
    main()
