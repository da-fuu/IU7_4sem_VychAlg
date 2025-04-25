from math import sqrt


def gauss(A, B):
    n = len(A)
    if n != len(B):
        raise ValueError("Матрица и вектор правой части должны иметь одинаковую длину")

    for i in range(n):
        for j in range(i + 1, n):
            coeff = -(A[j][i] / A[i][i])
            for k in range(i, n):
                A[j][k] += coeff * A[i][k]
            B[j] += coeff * B[i]

    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            coeff = -(A[j][i] / A[i][i])
            A[j][i] += coeff * A[i][i]
            B[j] += coeff * B[i]

    res = [B[i] / A[i][i] for i in range(n)]
    return res


def newton_system(j_func, funcs, start_root_approx, iter_limit, eps=1e-6):

    def f(x):
        return [func(*x) for func in funcs]

    xk = start_root_approx
    n = 1
    
    while True:
        dx = gauss(j_func(*xk), [-y for y in f(xk)])
        
        xnext = [xk[i] + dx[i] for i in range(len(xk))]

        if sqrt(sum([x**2 for x in dx])) < eps or n == iter_limit:
            return xnext, n
        
        xk = xnext
        n += 1


# Определение системы уравнений
def f1(x, y, z):
    return x ** 2 + y ** 2 + z ** 2 - 1


def f2(x, y, z):
    return 2 * x ** 2 + y ** 2 - 4 * z


def f3(x, y, z):
    return 3 * x ** 2 - 4 * y + z ** 2


def jacobian(x, y, z):
    return [  # J[i][j] = df_i(x_j) / dx_j
        [2 * x, 2 * y, 2 * z],
        [4 * x, 2 * y, -4],
        [6 * x, -4, 2 * z],
    ]


def solve_system():
    res, iters = newton_system(jacobian, [f1, f2, f3], [1, 1, 1], iter_limit=100)
    xres, yres, zres = res

    print("\n" + "=" * 50)
    print(f"{'Решение системы уравнений':^50}")
    print("-" * 50)
    print(f"{'Полученное решение':<30} | {'Значение':>15}")
    print("-" * 50)
    print(f"{'x':<30} | {xres:>15.6f}")
    print(f"{'y':<30} | {yres:>15.6f}")
    print(f"{'z':<30} | {zres:>15.6f}")
    print(f"{'Итераций':<30} | {iters:>15}")
    print("-" * 50)
