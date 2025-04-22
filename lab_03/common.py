
def solve_matrix_gauss(matrix):
    n = len(matrix)

    for i in range(n):
        for j in range(i + 1, n):
            if i == j:
                continue

            k = matrix[j][i] / matrix[i][i]

            for q in range(i, n + 1):
                matrix[j][q] -= k * matrix[i][q]

    result = [0 for _ in range(n)]

    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            matrix[i][n] -= result[j] * matrix[i][j]

        result[i] = matrix[i][n] / matrix[i][i]

    return result
