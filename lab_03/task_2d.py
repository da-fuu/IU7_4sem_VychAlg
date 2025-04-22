import numpy as np
import matplotlib.pyplot as plt
from random import random
from common import solve_matrix_gauss

SIZE_TABLE = 6


def generate_table_2d():
    table = []

    for i in range(SIZE_TABLE):
        table.append([random() * 5, random() * 5, random() * 5, random() * 5])

    return table


def solve_task_2d(table):
    xs = np.array([i[0] for i in table])
    ys = np.array([i[1] for i in table])
    zs = np.array([i[2] for i in table])
    weights = np.array([i[3] for i in table])

    sum_rho = np.sum(weights)
    sum_rho_x = np.sum(weights * xs)
    sum_rho_y = np.sum(weights * ys)
    sum_rho_z = np.sum(weights * zs)
    sum_rho_x2 = np.sum(weights * xs**2)
    sum_rho_y2 = np.sum(weights * ys**2)
    sum_rho_xy = np.sum(weights * (xs * ys))
    sum_rho_xz = np.sum(weights * (xs * zs))
    sum_rho_yz = np.sum(weights * (ys * zs))

    matrix = [
        [sum_rho_x2, sum_rho_xy, sum_rho_x, sum_rho_xz],
        [sum_rho_xy, sum_rho_y2, sum_rho_y, sum_rho_yz],
        [sum_rho_x, sum_rho_y, sum_rho, sum_rho_z]
    ]

    coefficients = solve_matrix_gauss(matrix)

    a, b, c = coefficients

    x_points = np.linspace(min(xs), max(xs), 100)
    y_points = np.linspace(min(ys), max(ys), 100)
    x_points, y_points = np.meshgrid(x_points, y_points)

    z_points = a * x_points + b * y_points + c

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(xs, ys, zs, c='r', marker='o', s=weights*10, label='Точки (размер = вес)')

    ax.plot_surface(x_points, y_points, z_points, alpha=0.5, color='b', label='Аппроксимирующая поверхность')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    plt.show()
