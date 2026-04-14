import os
import re
import numpy as np
import matplotlib.pyplot as plt


def read_data_file(filename):
    with open(filename, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    # число групп
    n_groups = int(lines[0])

    # количество точек в группах
    counts = list(map(int, lines[1].split()))
    if len(counts) != n_groups:
        raise ValueError("Количество групп не совпадает с числом элементов во 2-й строке")

    groups = []
    current_line = 2

    for g in range(n_groups):
        # пропускаем заголовок группы, если он есть
        group_name = f"Group_{g+1}"
        if current_line < len(lines) and lines[current_line].startswith("#"):
            header = lines[current_line]
            current_line += 1

            # пробуем вытащить имя элемента из заголовка
            # например: "#3 He--" -> "He"
            m = re.search(r'#\d+\s+([A-Za-z]+)', header)
            if m:
                group_name = m.group(1)

        points = []
        for _ in range(counts[g]):
            if current_line >= len(lines):
                raise ValueError("Файл 1 закончился раньше, чем ожидалось")

            parts = lines[current_line].split()
            current_line += 1

            if len(parts) < 4:
                raise ValueError(f"Некорректная строка данных: {lines[current_line-1]}")

            x = float(parts[2])
            y = float(parts[3])
            points.append((x, y))

        groups.append({
            "name": group_name,
            "points": points
        })

    return groups

def read_params_file(filename, n_groups):
    a = None
    b_dict = {}

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # a= ...
            m_a = re.match(r'^a=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)$', line)
            if m_a:
                a = float(m_a.group(1))
                continue

            # b( 1 )= ...
            m_b = re.match(
                r'^b\(\s*(\d+)\s*\)\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)$',
                line
            )
            if m_b:
                idx = int(m_b.group(1))
                val = float(m_b.group(2))
                b_dict[idx] = val
                continue

            # sigma ... пропускаем
            if line.lower().startswith("sigma"):
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

    for j, group in enumerate(groups):
        name = group["name"]
        points = group["points"]
        bj = b[j]

        x_data = np.array([p[0] for p in points])
        y_data = np.array([p[1] for p in points])

        # диапазон X для линии
        x_min = min(x_data.min(), 0) - 2
        x_max = max(x_data.max(), 0) + 2
        x_line = np.linspace(x_min, x_max, 200)
        y_line = a * x_line + bj

        plt.figure(figsize=(8, 5))
        plt.scatter(x_data, y_data, marker='s', color='black', label='Эксп. данные')
        plt.plot(x_line, y_line, color='black', linewidth=1.0,
                 label=f'y = {a:.2f}*x + {bj:.2f}')

        plt.title(name, fontsize=18, fontweight='bold')
        plt.xlabel('Qgg (МэВ)', fontsize=14, fontweight='bold', loc='right')
        plt.ylabel('Ln(σ)', fontsize=14, fontweight='bold', rotation=0, loc='top')

        plt.grid(False)
        plt.legend(frameon=False, loc='upper left')
        plt.tight_layout()

        out_path = os.path.join(output_dir, f"{j+1:02d}_{name}.png")
        plt.savefig(out_path, dpi=150)
        plt.close()

        print(f"Сохранен график: {out_path}")


def main():
    data_file = "../input/OTa_qgg.dat"
    params_file = "../output/OTa_qgg_out.dat"

    groups = read_data_file(data_file)
    a, b = read_params_file(params_file, len(groups))

    print(f"a = {a}")
    for i, bj in enumerate(b, start=1):
        print(f"b{i} = {bj}")

    plot_groups(groups, a, b, output_dir="plots")


if __name__ == "__main__":
    main()
