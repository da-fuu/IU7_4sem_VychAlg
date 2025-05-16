import numpy as np

def read_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    x_line = lines[0].strip().split('\t')
    x_values = [float(val) for val in x_line[1:] if val]

    y_values = []
    z_matrix = []
    for line in lines[1:]:
        parts = line.strip().split('\t')
        y_values.append(float(parts[0]))
        z_row = [float(val) for val in parts[1:] if val]
        z_matrix.append(z_row)

    return np.array(x_values), np.array(y_values), np.array(z_matrix)


def bilinear_interpolation(x, y, x_grid, y_grid, z_matrix):
    # Находим ближайшие индексы
    i = np.searchsorted(x_grid, x) - 1
    j = np.searchsorted(y_grid, y) - 1

    # Ограничиваем индексы в пределах массива
    i = max(0, min(i, len(x_grid) - 2))
    j = max(0, min(j, len(y_grid) - 2))

    # Координаты соседних точек
    x0, x1 = x_grid[i], x_grid[i + 1]
    y0, y1 = y_grid[j], y_grid[j + 1]

    # Значения функции в соседних точках
    f00 = z_matrix[j][i]
    f01 = z_matrix[j][i + 1]
    f10 = z_matrix[j + 1][i]
    f11 = z_matrix[j + 1][i + 1]

    # Вычисляем веса
    dx = x1 - x0 if x1 != x0 else 1.0
    dy = y1 - y0 if y1 != y0 else 1.0
    wx = (x - x0) / dx
    wy = (y - y0) / dy

    # Билинейная интерполяция
    return (f00 * (1 - wx) * (1 - wy) +
            f01 * wx * (1 - wy) +
            f10 * (1 - wx) * wy +
            f11 * wx * wy)


def double_integrate(x_grid, y_grid, z_matrix, phi, psi, a, b, n_x=100, n_y=100):
    x_nodes = np.linspace(a, b, n_x)
    integral = 0.0

    for x in x_nodes:
        y_a = phi(x)
        y_b = psi(x)
        if y_a >= y_b:
            continue

        y_nodes = np.linspace(y_a, y_b, n_y)
        f_values = np.array([bilinear_interpolation(x, y, x_grid, y_grid, z_matrix) for y in y_nodes])

        # Метод трапеций по y
        h_y = (y_b - y_a) / (n_y - 1)
        integral_y = h_y * (0.5 * f_values[0] + 0.5 * f_values[-1] + np.sum(f_values[1:-1]))
        integral += integral_y

    # Умножаем на шаг по x
    h_x = (b - a) / (n_x - 1)
    integral *= h_x

    return integral


def task_2():
    file_path = "data.txt"
    x_values, y_values, z_matrix = read_data(file_path)

    alpha = 1.0
    beta = 4.0
    a = 0.0
    b = 2.0

    # alpha, beta, a, b = map(float, input('Введите Альфа, Бета, a, b: ').split())

    def phi(x):
        return alpha * x ** 2

    def psi(x):
        return beta * x ** 2

    integral = double_integrate(x_values, y_values, z_matrix, phi, psi, a, b)
    print(f"Двойной интеграл: {integral}")
    45