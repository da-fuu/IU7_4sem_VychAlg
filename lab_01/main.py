from math import factorial
from copy import deepcopy

EPS = 1e-10


class Interpolator:
    def __init__(self, filename):
        self.table = []
        with open(filename) as f:
            f.readline()
            for s in f:
                self.table.append([float(d) for d in s.split()])
        self.table.sort(key=lambda line: line[0])

    def reverse_table(self):
        new_table = [[] for _ in range(len(self.table))]
        for i in range(len(self.table)):
            new_table[i].append(self.table[i][1])
            new_table[i].append(self.table[i][0])
            if len(self.table[0]) == 4:
                new_table[i].append(1 / self.table[i][2])
                new_table[i].append(-self.table[i][3] / (self.table[i][2] ** 3))
        self.table = new_table
        self.table.sort(key=lambda line: line[0])

    def find_start_index(self, x: float, points: int):
        after = -1
        while after + 1 < len(self.table) and x > self.table[after + 1][0]:
            after += 1
        if after < points - 1:
            return 0
        if after >= len(self.table) - points:
            return len(self.table) - points
        return after - points // 2 + 1

    def find_diff_derivative(self, index, order):
        return self.table[index][order + 1] / factorial(order)

    def create_diffs(self, indexes):
        diffs = [[] for _ in range(len(indexes))]
        for ind_el in indexes:
            diffs[0].append(self.table[ind_el][1])
        xs = list(zip(*self.table))[0]
        for order in range(1, len(indexes)):
            for ind_el in range(len(indexes) - order):
                divisor = xs[indexes[ind_el]] - xs[indexes[ind_el + order]]
                if abs(divisor) < EPS:
                    diffs[order].append(self.find_diff_derivative(indexes[ind_el], order))
                else:
                    diffs[order].append((diffs[order-1][ind_el] - diffs[order-1][ind_el + 1]) / divisor)
        return [diff[0] for diff in diffs]

    def find_polynom_val(self, x: float, indexes, diffs):
        ans = 0
        for i in range(len(indexes)):
            mul = 1
            for j in range(i):
                mul *= x - self.table[indexes[j]][0]
            ans += diffs[i] * mul
        return ans

    def interpolate_newton(self, x: float, n: int):
        if n < 0:
            print('Введена отрицательная степень полинома Ньютона!')
            return
        if n >= len(self.table):
            print(f'Введена степень полинома Ньютона больше {len(self.table) - 1}!')
            return
        start_index = self.find_start_index(x, n + 1)
        indexes = list(range(start_index, start_index + n + 1))
        diffs = self.create_diffs(indexes)
        y = self.find_polynom_val(x, indexes, diffs)
        return y

    def interpolate_ermit(self, x: float, points: int):
        if points < 1:
            print('Введено количество узлов полинома Эрмита меньше 1!')
            return
        if points > len(self.table):
            print(f'Введено количество узлов полинома Эрмита больше {len(self.table)}!')
            return
        start_index = self.find_start_index(x, points)
        indexes = [i // 3 for i in range(start_index * 3, (start_index + points) * 3)]
        diffs = self.create_diffs(indexes)
        y = self.find_polynom_val(x, indexes, diffs)
        return y

    def find_next_index(self, x: float) -> int:
        for i in range(len(self.table)):
            if self.table[i][0] > x:
                return i - 1

    def find_run_k(self) -> list[list[float]]:
        n: int = len(self.table) - 1
        y = list(zip(*self.table))[1]
        z: list[float] = [100000 for _ in range(n)]
        m: list[float] = [100000 for _ in range(n)]
        h = 0.2

        z[0] = 0
        m[0] = 0

        for i in range(2, n + 1):
            z[i - 1] = -h / (h * z[i - 2] + 4 * h)
            f = 3 * ((y[i] - y[i - 1]) / h - (y[i - 1] - y[i - 2]) / h)
            m[i - 1] = (f - h * m[i - 2]) / (h * z[i - 2] + 4 * h)
        return list(zip(z, m))

    def calc_splines(self):
        run_k: list[list[float]] = self.find_run_k()

        n: int = len(self.table) - 1
        y = list(zip(*self.table))[1]
        a: list[float] = [10000 for _ in range(n)]
        b: list[float] = [10000 for _ in range(n)]
        c: list[float] = [10000 for _ in range(n + 1)]
        d: list[float] = [10000 for _ in range(n)]
        h = 0.2

        c[0] = 0
        c[n] = 0
        for i in range(n, 1, -1):
            c[i - 1] = run_k[i - 2][0] * c[i] + run_k[i - 2][1]

        for i in range(1, n):
            b[i - 1] = (y[i] - y[i - 1]) / h - h * (c[i] - 2 * c[i - 1]) / 3
        b[n - 1] = (y[n] - y[n - 1]) / h - h * 2 * c[n - 1] / 3

        for i in range(1, n + 1):
            a[i - 1] = y[i - 1]

        d[n - 1] = -c[n - 1] / 3 / h
        for i in range(1, n):
            d[i - 1] = (c[i] - c[i - 1]) / 3 / h
        return list(zip(a, b, c[:n], d))

    def interpolate_spline(self, x: float) -> float:
        splines: list[list[float]] = self.calc_splines()
        ind: int = self.find_next_index(x)
        k: list[float] = list(splines[ind])
        k[2] *= splines[ind+1][2]
        xi: float = self.table[ind][0]
        y: float = k[0] + k[1] * (x - xi) + k[2] * (x - xi) ** 2 + k[3] * (x - xi) ** 3
        return y


def main():
    task = 1
    # task = int(input('Введите номер решаемой задачи (1 - интерполяция в точке обоими полиномами,'
    #                  ' 2 - сравнение полиномов, 3 - поиск корня функции, 4 - решение системы уравнений):\n'))
    if not 1 <= task <= 4:
        print('Введен неверный номер задания!')
        return

    if task == 1:
        interpolator = Interpolator('data1.txt')

        import matplotlib.pyplot as plt
        import numpy as np

        xs = np.linspace(0.1, 4, 600)
        ys = []
        ys2 = []
        for x in xs:
            ys.append(interpolator.interpolate_spline(x))
            ys2.append(interpolator.interpolate_newton(x, 5))

        plt.plot(xs, ys)
        plt.plot(xs, ys2)
        plt.show()

        n = int(input('Введите степень полинома Ньютона (от 0 до количества точек в файле - 1): '))
        points = int(input('Введите количество узлов для полинома Эрмита (от 1 до количества точек в файле): '))
        x = float(input('Введите значение х для интерполяции: '))
        y_newton = interpolator.interpolate_newton(x, n)
        y_ermit = interpolator.interpolate_ermit(x, points)
        if y_newton is None or y_ermit is None:
            return
        print('Полученное значение полиномом Ньютона: {:g}'.format(y_newton))
        print('Полученное значение полиномом Эрмита: {:g}'.format(y_ermit))

    elif task == 2:
        interpolator = Interpolator('data1.txt')
        x = float(input('Введите значение х для сравнения полиномов: '))
        results = [dict() for _ in range(1, 6)]
        for points in range(1, 6):
            results[points - 1]['n'] = points * 3 - 1
            results[points - 1]['newton'] = interpolator.interpolate_newton(x, points * 3 - 1)
            results[points - 1]['ermit'] = interpolator.interpolate_ermit(x, points)
        print('{:7s} | {:8s} | {:8s} |'.format('Степень', 'Ньютон', 'Эрмит'))
        for res in results:
            print('{:7d} | {:8g} | {:8g} |'.format(res['n'], res['newton'], res['ermit']))

    elif task == 3:
        interpolator = Interpolator('data1.txt')
        interpolator.table = interpolator.table[:6]
        n = int(input('Введите степень полинома Ньютона (от 0 до 5): '))
        points = int(input('Введите количество узлов для полинома Эрмита (от 1 до 6): '))
        interpolator.reverse_table()
        root_newton = interpolator.interpolate_newton(0, n)
        root_ermit = interpolator.interpolate_ermit(0, points)
        if root_newton is None or root_ermit is None:
            return
        print('Корень уравнения рассчитанный полиномом Ньютона: {:g}'.format(root_newton))
        print('Корень уравнения рассчитанный полиномом Эрмита: {:g}'.format(root_ermit))

    elif task == 4:
        n = int(input('Введите степень полинома Ньютона (от 0 до 5): '))
        interpolator = Interpolator('data2_1.txt')
        interpolator.reverse_table()
        temp_interpolator = Interpolator('data2_2.txt')

        new_table = deepcopy(interpolator.table)
        for i in range(len(new_table)):
            inter_y = temp_interpolator.interpolate_newton(new_table[i][0], n)
            if inter_y is None:
                return
            new_table[i][1] -= inter_y
        interpolator.table = new_table

        interpolator.reverse_table()
        interpolator.table = interpolator.table[:6]
        root_x = interpolator.interpolate_newton(0, n)
        root_y = temp_interpolator.interpolate_newton(root_x, n)
        print('Корень системы уравнений рассчитанный полиномом Ньютона: {{{:g}, {:g}}}'.format(root_x, root_y))


if __name__ == '__main__':
    main()
