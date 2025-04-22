from random import random
from copy import deepcopy
import matplotlib.pyplot as plt
from common import solve_matrix_gauss

SIZE_TABLE = 6
EPS = 0.01


def generate_table():
    table = []
    x = random() * 5

    for i in range(SIZE_TABLE):
        x += random() + 0.5
        table.append([x, random() * 5, 1])

    return table


def table_changed(table):
    for i in table:
        if i[2] != 1:
            return True

    return False


def get_base_table(table):
    base_table = deepcopy(table)

    for i in range(len(base_table)):
        base_table[i][2] = 1

    return base_table


def make_sle_matrix(table, n):
    size = len(table)
    matrix = [[0 for _ in range(n + 2)] for _ in range(n + 1)]

    for i in range(n + 1):
        for j in range(n + 1):
            matrix[i][j] = 0.0
            matrix[i][n + 1] = 0.0

            for k in range(size):
                weight = table[k][2]
                x = table[k][0]
                y = table[k][1]

                matrix[i][j] += weight * pow(x, (i + j))    # Сумма ρ_k * x_k^(i+j)
                matrix[i][n + 1] += weight * y * pow(x, i)  # Сумма ρ_k * y_k * x_k^i

    return matrix


def find_graph_dots(table, n):
    matrix = make_sle_matrix(table, n)
    result = solve_matrix_gauss(matrix)

    x_arr = []
    y_arr = []
    k = table[0][0] - EPS

    while k <= table[len(table) - 1][0] + EPS:
        y = 0
        for j in range(0, n + 1):
            y += result[j] * pow(k, j)  # a * k ** j

        x_arr.append(k)
        y_arr.append(y)

        k += EPS

    return x_arr, y_arr


def plot_graphs(table, n, type_graph, type_dots):
    for i in sorted({1, 2, n}):
        x_arr, y_arr = find_graph_dots(table, i)
        plt.plot(x_arr, y_arr, type_graph, label="{:s}\nn = {:d}".format(type_dots, i))



def solve_task_1d(table):
    try:
        n = int(input("\nВведите степень аппроксимирующего полинома: "))
    except ValueError:
        print("\nОшибка: некорректно введена степень полинома!")
        return

    if n < 0:
        print("\nОшибка: некорректно введена степень полинома!")
        return

    if table_changed(table):
        base_table = get_base_table(table)
        type_dots = "Diff weights"
        type_graph = "-."

        plot_graphs(base_table, n, "-", "Equal weights")
    else:
        type_dots = "Equal weights"
        type_graph = "-"

    plot_graphs(table, n, type_graph, type_dots)

    x_arr = [i[0] for i in table]
    y_arr = [i[1] for i in table]

    plt.plot(x_arr, y_arr, 'o')

    plt.legend()
    plt.grid()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()
