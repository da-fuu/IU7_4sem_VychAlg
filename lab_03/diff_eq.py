import numpy as np
import matplotlib.pyplot as plt
from common import solve_matrix_gauss

# Базисные функции
def u0(x):
    return 1 - x  # Удовлетворяет краевым условиям y(0) = 1, y(1) = 0

def u1(x):
    return x * (1 - x)

def u2(x):
    return x ** 2 * (1 - x)

def u3(x):
    return x ** 3 * (1 - x)

u = [u0, u1, u2, u3]


# Производные базисных функций
def du0(x):
    return -1

def du1(x):
    return 1 - 2 * x

def du2(x):
    return 2 * x - 3 * x ** 2

def du3(x):
    return 3 * x ** 2 - 4 * x ** 3

du = [du0, du1, du2, du3]


# Вторая производная базисных функций
def d2u0(x):
    return 0

def d2u1(x):
    return -2

def d2u2(x):
    return 2 - 6 * x

def d2u3(x):
    return 6 * x - 12 * x ** 2

d2u = [d2u0, d2u1, d2u2, d2u3]


# Правая часть
def f(x):
    return 2 * x


def trapezoidal_integral(func, a, b, n=1000):
    h = (b - a) / n
    integral = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        integral += func(a + i * h)
    integral *= h
    return integral


# Формируем систему уравнений для m
def solve_for_m(m):
    A = [[0 for _ in range(m)] for _ in range(m)]
    B = [0 for _ in range(m)]

    for i in range(1, m + 1):
        for j in range(1, m + 1):
            def integral_left(x):
                return (d2u[i](x) + x * du[i](x) + u[i](x)) * u[j](x)

            A[i - 1][j - 1] = trapezoidal_integral(integral_left, 0, 1)

        def integral_right(x):
            return (f(x) - (d2u0(x) + x * du0(x) + u0(x))) * u[i](x)

        B[i - 1] = trapezoidal_integral(integral_right, 0, 1)

    matrix = [row + [B[idx]] for idx, row in enumerate(A)]

    return solve_matrix_gauss(matrix)


def solve_diff_eq():
    C_m2 = solve_for_m(2)
    print("Коэффициенты для m=2:", C_m2)

    C_m3 = solve_for_m(3)
    print("Коэффициенты для m=3:", C_m3)

    x = np.linspace(0, 1, 100)
    y_m2 = u0(x) + C_m2[0] * u1(x) + C_m2[1] * u2(x)
    plt.plot(x, y_m2, label='m=2')
    y_m3 = u0(x) + C_m3[0] * u1(x) + C_m3[1] * u2(x) + C_m3[2] * u3(x)
    plt.plot(x, y_m3, label='m=3')

    plt.legend()
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Сравнение решений для m=2 и m=3')
    plt.show()
