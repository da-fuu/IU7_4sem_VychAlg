length = 5


def read_data(filename: str) -> list[list[list[float]]]:
    table: list[list[list[float]]] = []
    with open(filename) as f:
        for z in range(length):
            f.readline()
            f.readline()
            table.append([list(map(float, f.readline().split()[1:])) for _ in range(length)])
            f.readline()
    return table


class Interpolator:
    def __init__(self, table: list[list[float]]):
        self.table: list[list[float]] = table

    def find_start_index(self, x: float, points: int) -> int:
        after: int = -1
        while after + 1 < len(self.table) and x > self.table[after + 1][0]:
            after += 1
        if after < points - 1:
            return 0
        if after >= len(self.table) - points:
            return len(self.table) - points
        return after - points // 2 + 1

    def create_diffs(self, indexes: list[int]) -> list[float]:
        diffs: list[list[float]] = [[] for _ in range(len(indexes))]
        for ind_el in indexes:
            diffs[0].append(self.table[ind_el][1])
        xs: list[float] = list(zip(*self.table))[0]
        for order in range(1, len(indexes)):
            for ind_el in range(len(indexes) - order):
                divisor: float = xs[indexes[ind_el]] - xs[indexes[ind_el + order]]
                diffs[order].append((diffs[order-1][ind_el] - diffs[order-1][ind_el + 1]) / divisor)
        return [diff[0] for diff in diffs]

    def find_polynom_val(self, x: float, indexes: list[int], diffs: list[float]) -> float:
        ans: float = 0
        for i in range(len(indexes)):
            mul: float = 1
            for j in range(i):
                mul *= x - self.table[indexes[j]][0]
            ans += diffs[i] * mul
        return ans

    def interpolate_newton(self, x: float, n: int) -> float:
        start_index: int = self.find_start_index(x, n + 1)
        indexes: list[int] = list(range(start_index, start_index + n + 1))
        diffs: list[float] = self.create_diffs(indexes)
        y: float = self.find_polynom_val(x, indexes, diffs)
        return y

    def find_prev_index(self, x: float) -> int:
        for i in range(len(self.table)):
            if self.table[i][0] > x:
                return i - 1

    def find_run_k(self) -> list[list[float]]:
        n: int = len(self.table) - 1
        y = list(zip(*self.table))[1]
        z: list[float] = [100000 for _ in range(n)]
        m: list[float] = [100000 for _ in range(n)]
        h = self.table[1][0] - self.table[0][0]

        z[0] = 0
        m[0] = 0

        for i in range(2, n + 1):
            z[i - 1] = - h / (h * z[i - 2] + 4 * h)
            f = 3 * ((y[i] - y[i - 1]) / h - (y[i - 1] - y[i - 2]) / h)
            m[i - 1] = (f - h * m[i - 2]) / (h * z[i - 2] + 4 * h)
# 3 4 5
        return list(zip(z, m))

    def calc_splines(self):
        run_k: list[list[float]] = self.find_run_k()

        n: int = len(self.table) - 1
        y = list(zip(*self.table))[1]
        a: list[float] = [10000 for _ in range(n)]
        b: list[float] = [10000 for _ in range(n)]
        c: list[float] = [10000 for _ in range(n + 1)]
        d: list[float] = [10000 for _ in range(n)]
        h = self.table[1][0] - self.table[0][0]

        c[0] = 0.0
        c[n] = run_k[-1][1]
        for i in range(n, 1, -1):
            c[i - 1] = run_k[i - 2][0] * c[i] + run_k[i - 2][1]
        c = c[1:]

        for i in range(1, n):
            b[i - 1] = (y[i] - y[i - 1]) / h - h * (c[i] + 2 * c[i - 1]) / 3
        b[n - 1] = (y[n] - y[n - 1]) / h - h * 2 * c[n - 1] / 3

        for i in range(1, n + 1):
            a[i - 1] = y[i - 1]

        d[n - 1] = -c[n - 1] / 3 / h
        for i in range(1, n):
            d[i - 1] = (c[i] - c[i - 1]) / 3 / h
        return list(zip(a, b, c, d))

    def interpolate_spline(self, x: float) -> float:
        splines: list[list[float]] = self.calc_splines()
        ind: int = self.find_prev_index(x)
        k: list[float] = splines[ind]
        xi: float = self.table[ind][0]
        y: float = k[0] + k[1] * (x - xi) + k[2] * (x - xi) ** 2 + k[3] * (x - xi) ** 3
        return y


def get_data_0d(data_1d: list[float], x: float, n: int) -> float:
    table: list[list[float]] = [[float(i), data_1d[i]] for i in range(length)]
    inter: Interpolator = Interpolator(table)
    if n >= 0:
        return inter.interpolate_newton(x, n)
    else:
        return inter.interpolate_spline(x)


def get_data_1d(data_2d, y: float, n) -> list[float]:
    data_1d: list[float] = []
    for x in range(length):
        table = [data_2d[i][x] for i in range(length)]
        data_1d.append(get_data_0d(table, y, n))
    return data_1d


def get_data_2d(data_3d, z: float, n) -> list[list[float]]:
    data_2d: list[list[float]] = []
    for y in range(length):
        table = [data_3d[i][y] for i in range(length)]
        data_2d.append(get_data_1d(table, z, n))
    return data_2d


def main():
    data_3d: list[list[list[float]]] = read_data('data.txt')
    x, y, z = map(float, input('Введите координаты точки x, y, z:\n').split())
    nx, ny, nz = map(int, input('Введите степени полинома Ньютона по осям x, y, z или -1 для интерполяции сплайнами '
                                '(целое число от -1 до 4):\n').split())
    data_2d: list[list[float]] = get_data_2d(data_3d, z, nz)
    data_1d: list[float] = get_data_1d(data_2d, y, ny)
    data_0d: float = get_data_0d(data_1d, x, nx)
    print('Интерполированное значение: {:g}'.format(data_0d))


if __name__ == '__main__':
    main()
