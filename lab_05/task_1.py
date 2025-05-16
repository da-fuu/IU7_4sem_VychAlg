import numpy as np
from math import fabs

def roots_legendre(n):
    if n < 1:
        return 

    if n == 1:
        return np.array([0.0]), np.array([2.0])

    beta = [k / np.sqrt(4*k**2 - 1) for k in range(1, n)]
    J = np.diag(beta, 1) + np.diag(beta, -1)

    eigenvalues, eigenvectors = np.linalg.eigh(J)

    nodes = eigenvalues
    weights = 2 * (eigenvectors[0, :])**2 

    idx = np.argsort(nodes)
    nodes = nodes[idx]
    weights = weights[idx]

    return nodes, weights

def integrate_trapezoid(f, a, b, n=100):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    integral = h * (0.5 * y[0] + 0.5 * y[-1] + np.sum(y[1:-1]))
    return integral


def integrate_simpson(f, a, b, n=100):
    if n % 2 != 0:
        n += 1  # Делаем n четным
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    integral = h / 3 * (y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-2:2]))
    return integral


def integrate_gauss(f, a, b, n=100):
    # Получаем узлы и веса для полинома Лежандра
    nodes, weights = roots_legendre(n)
    # Преобразуем узлы из интервала [-1, 1] в [a, b]
    transformed_nodes = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    integral = 0.5 * (b - a) * np.sum(weights * f(transformed_nodes))
    return integral


def func_k1(x):
    return abs(x)

def func_k2(x):
    return x ** 2

def task_1():   
    # Аналитические решения для сравнения
    exact_k1 = 1.0
    exact_k2 = 2 / 3

    a, b = -1, 1  # Границы интегрирования
    n = 3  # Количество узлов

    print("\nРезультаты для k=1 (функция |x|):")
    print(f"Точное значение: {exact_k1:.6f}")
    print(
        f"Метод трапеций: {integrate_trapezoid(func_k1, a, b, n):.6f}, ошибка: {fabs(integrate_trapezoid(func_k1, a, b, n) - exact_k1):.6f}")
    print(
        f"Метод Симпсона: {integrate_simpson(func_k1, a, b, n):.6f}, ошибка: {fabs(integrate_simpson(func_k1, a, b, n) - exact_k1):.6f}")
    print(
        f"Метод Гаусса: {integrate_gauss(func_k1, a, b, n):.6f}, ошибка: {fabs(integrate_gauss(func_k1, a, b, n) - exact_k1):.6f}")

    print("\nРезультаты для k=2 (функция x**2):")
    print(f"Точное значение: {exact_k2:.6f}")
    print(
        f"Метод трапеций: {integrate_trapezoid(func_k2, a, b, n):.6f}, ошибка: {fabs(integrate_trapezoid(func_k2, a, b, n) - exact_k2):.6f}")
    print(
        f"Метод Симпсона: {integrate_simpson(func_k2, a, b, n):.6f}, ошибка: {fabs(integrate_simpson(func_k2, a, b, n) - exact_k2):.6f}")
    print(
        f"Метод Гаусса: {integrate_gauss(func_k2, a, b, n):.6f}, ошибка: {fabs(integrate_gauss(func_k2, a, b, n) - exact_k2):.6f}\n")
