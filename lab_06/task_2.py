import numpy as np
import matplotlib.pyplot as plt


def solve_boundary_problem(alpha, beta, gamma, n=100):
    # Параметры сетки
    a = 0.0
    b = 1.0
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)

    # Коэффициенты разностной схемы
    A = np.zeros(n + 1)
    B = np.zeros(n + 1)
    C = np.zeros(n + 1)
    F = np.zeros(n + 1)

    # Заполнение коэффициентов для внутренних узлов
    for i in range(1, n):
        xi = x[i]
        A[i] = 1 / h ** 2 - xi ** 2 / h
        B[i] = -2 / h ** 2 + 4
        C[i] = 1 / h ** 2 + xi ** 2 / h
        F[i] = 2 * xi + np.exp(-xi)

    # Левое граничное условие.
    # Используем аппроксимацию второго порядка точности
    A[0] = 0
    B[0] = -3 / (2 * h)
    C[0] = 4 / (2 * h)
    F[0] = -1 / (2 * h) + alpha

    # Правое граничное условие.
    # Используем аппроксимацию второго порядка точности
    A[n] = -4 / (2 * h)
    B[n] = 3 / (2 * h) - beta
    C[n] = 0
    F[n] = 1 / (2 * h) - gamma

    u = thomas_algorithm(A, B, C, F)

    return x, u


def thomas_algorithm(a, b, c, d):
    n = len(d)

    # Прямой ход прогонки
    cp = np.zeros(n)
    dp = np.zeros(n)
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]

    for i in range(1, n):
        denom = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / denom if i < n - 1 else 0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom

    # Обратный ход прогонки
    x = np.zeros(n)
    x[-1] = dp[-1]

    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]

    return x


def task2():
    # Параметры задачи
    alpha = 1.0
    beta = 0.5
    gamma = 0.2

    x, u = solve_boundary_problem(alpha, beta, gamma)

    for i in range(len(x)):
        print(f'x = {x[i]:.2f}, u = {u[i]:.6f}')

    plt.plot(x, u)
    plt.show()
