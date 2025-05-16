import numpy as np
import math


def sympson_integral(x, f):
    def symfunc(a: float, b: float) -> float:
        return ((b - a) / 6) * (f(a) + 4 * f((a + b) / 2) + f(b))

    return sum([symfunc(x[i], x[i + 1]) for i in range(len(x) - 1)])


def half_division(func, target, start, end, iter_limit=10, eps=1e-6):
    for i in range(1, iter_limit + 1):
        center = (start + end) / 2

        if abs(func(center) - target) < eps or i == iter_limit:
            return center, i

        if (func(start) - target) * (func(center) - target) < 0:
            end = center
        else:
            start = center


# Функция для вычисления функции Лапласа (нормального распределения)
def laplas(x, integral_count=10):

    def under_integral_func(t):
        return np.exp(-(t ** 2) / 2)

    linspace_to_x = list(np.linspace(0, x, integral_count))
    return 2 / np.sqrt(2 * math.pi) * sympson_integral(linspace_to_x, under_integral_func)


# Функция для нахождения аргументов функции Лапласа
def find_laplace_argument():
    x = np.linspace(-5, 5, 100)

    target_y = [0, 0.14, -0.7, 0.90]
    found_x = [half_division(laplas, y_val, min(x), max(x), iter_limit=30)[0] for y_val in target_y]

    print("\n" + "=" * 40)
    print(f"{'Аргументы функции Лапласа':^40}")
    print("=" * 40)
    print(f"{'Заданное y':<20} | {'Найденное x':>15}")
    print("-" * 40)
    for y_val, x_val in zip(target_y, found_x):
        print(f"{y_val:<20.6f} | {x_val:>15.6f}")
    print("=" * 40)

    # plt.figure(figsize=(10, 6))
    # plt.plot(x, y, label='Функция Лапласа')
    # plt.scatter(found_x, target_y, color='red', label='Найденные точки')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.title('Функция Лапласа и найденные точки')
    # plt.legend()
    # plt.grid(True)
    # plt.show()
