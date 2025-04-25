import numpy as np
import matplotlib.pyplot as plt
from system_solve import newton_system

# Функции для решения краевой задачи
x0, y0 = 0, 1  # Граничное условие при x=0
x1, y1 = 1, 3  # Граничное условие при x=1
N = 100
step = (x1 - x0) / N


def jacobian_diff(*y):
    # Матрица Якоби для дискретизированного дифференциального уравнения
    n = len(y)
    res = []

    res.append([1] + [0] * (n - 1))  # Граничное условие при x=0

    for i in range(1, n - 1):
        row = [0] * (i - 1) + [1 / step ** 2] + [-2 / step ** 2 - 3 * y[i] ** 2] + [1 / step ** 2] + [0] * (n - i - 2)
        res.append(row)

    res.append([0] * (n - 1) + [1])  # Граничное условие при x=1
    return res


def f(n, x):
    # Определение функций для дискретизированного дифференциального уравнения
    if n == 0:
        def resf(*y: list[float | int]) -> float:
            return y[0] - y0
    elif n == N:
        def resf(*y: list[float | int]) -> float:
            return y[n] - y1
    else:
        def resf(*y: list[float | int]) -> float:
            return (y[n - 1] + -2 * y[n] + y[n + 1]) / step ** 2 - y[n] ** 3 - x[n] ** 2
    return resf


def starty(x):
    return 2 * x + 1


def solve_boundary_problem():
    x = np.linspace(x0, x1, N + 1)
    y = [starty(xp) for xp in x]

    funcs = [f(n, x) for n in range(N + 1)]

    # Решение методом Ньютона
    res, iters = newton_system(jacobian_diff, funcs, y, iter_limit=30, eps=1e-15)

    # Форматированный вывод результатов
    print("\n" + "-" * 40)
    print(f"{'Решение краевой задачи':^40}")
    print("-" * 40)
    print(f"{'Количество итераций':<20} | {iters:>15}")
    print("-" * 40)

    # Построение графика (если включено)
    plt.figure(figsize=(10, 6))
    plt.plot(x, res, label='Решение краевой задачи')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Численное решение краевой задачи y\'\' - y^3 = x^2')
    plt.legend()
    plt.grid(True)
    plt.show()
