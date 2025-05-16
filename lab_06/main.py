import numpy as np
import matplotlib.pyplot as plt  # matplotlib для графиков


# --- Метод прогонки (алгоритм Томаса) ---
def thomas_algorithm(a_coeff, b_coeff, c_coeff, d_rhs):
    n = len(d_rhs)
    c_prime = np.zeros(n)  # Временный массив для коэффициентов c
    d_prime = np.zeros(n)  # Временный массив для коэффициентов d
    x_solution = np.zeros(n)  # Вектор решения

    if n == 0:
        return x_solution
    if n == 1:
        if abs(b_coeff[0]) < 1e-14:  # Проверка деления на ноль
            raise ValueError("Метод прогонки не удался: b_coeff[0] равен нулю для одного уравнения.")
        x_solution[0] = d_rhs[0] / b_coeff[0]
        return x_solution

    # Прямой ход (вычисление прогоночных коэффициентов)
    if abs(b_coeff[0]) < 1e-14:  # Проверка деления на ноль
        raise ValueError("Метод прогонки не удался: b_coeff[0] равен нулю.")
    c_prime[0] = c_coeff[0] / b_coeff[0]
    d_prime[0] = d_rhs[0] / b_coeff[0]

    for i in range(1, n):
        denominator = b_coeff[i] - a_coeff[i] * c_prime[i - 1]
        if abs(denominator) < 1e-14:  # Знаменатель равен нулю или слишком мал
            raise ValueError(
                f"Метод прогонки не удался: нулевой знаменатель на шаге i={i} (b_coeff[{i}] - a_coeff[{i}]*c_prime[{i - 1}])")
        if i < n - 1:  # c_coeff[n-1] - последний элемент c_coeff, не используется для c_prime[n-1]
            c_prime[i] = c_coeff[i] / denominator
        d_prime[i] = (d_rhs[i] - a_coeff[i] * d_prime[i - 1]) / denominator

    # Обратный ход (вычисление решения)
    x_solution[n - 1] = d_prime[n - 1]
    for i in range(n - 2, -1, -1):
        x_solution[i] = d_prime[i] - c_prime[i] * x_solution[i + 1]

    return x_solution


# --- Задание 1: Численное дифференцирование ---
def task1_numerical_differentiation(x_data, y_data):
    """
    Выполняет численное дифференцирование для таблично заданной функции.
    Печатает результаты напрямую.
    """
    num_points = len(x_data)
    if num_points < 3:  # Нужно как минимум 3 точки для некоторых методов (например, центральная разность или Рунге с 2h)
        raise ValueError("Необходимо как минимум 3 точки данных для всех методов дифференцирования.")

    # Проверка равномерности шага
    h_steps = np.diff(x_data)
    if not np.allclose(h_steps, h_steps[0]):
        raise ValueError("Данные x_data должны иметь равномерный шаг для этих формул.")
    h = h_steps[0]

    # Инициализация массивов результатов с np.nan
    y_fwd_diff = np.full(num_points, np.nan)  # (1) Односторонняя разностная производная (вперед)
    y_central_diff = np.full(num_points, np.nan)  # (2) Центральная разностная производная
    y_runge = np.full(num_points, np.nan)  # (3) 2-я формула Рунге
    y_rectifying_vars = np.full(num_points, np.nan)  # (4) Выравнивающие переменные
    y_second_deriv = np.full(num_points, np.nan)  # (5) Вторая разностная производная

    # 1. Односторонняя разностная производная (вперед) (y_fwd_diff (1))
    # Определена для i = 0, ..., num_points-2
    for i in range(num_points - 1):
        y_fwd_diff[i] = (y_data[i + 1] - y_data[i]) / h

    # 2. Центральная разностная производная (y_central_diff (2))
    # Определена для i = 1, ..., num_points-2
    for i in range(1, num_points - 1):
        y_central_diff[i] = (y_data[i + 1] - y_data[i - 1]) / (2 * h)

    # 3. 2-я формула Рунге (используя односторонние разности вперед) (y_runge (3))
    # p=1 (порядок точности односторонней разности), m=2 (множитель шага)
    # Определена для i = 0, ..., num_points-3 (требует y[i], y[i+1], y[i+2])
    if num_points >= 3:
        for i in range(num_points - 2):
            d_h = (y_data[i + 1] - y_data[i]) / h
            d_2h = (y_data[i + 2] - y_data[i]) / (2 * h)
            y_runge[i] = d_h + (d_h - d_2h)  # Знаменатель (m^p-1) равен 1

    # 4. Выравнивающие переменные (y_rectifying_vars (4))
    # y = a0*x / (a1 + a2*x)  =>  1/y = a1/(a0*x) + a2/a0
    # eta = 1/y, xi = 1/x. Тогда eta = K*xi + M. d(eta)/d(xi) = K
    # y'_x = K * (y/x)^2
    # Определена для i = 0, ..., num_points-2

    safe_x_data = np.where(np.abs(x_data) < 1e-12, np.nan, x_data)
    safe_y_data = np.where(np.abs(y_data) < 1e-12, np.nan, y_data)

    xi_values = 1.0 / safe_x_data
    eta_values = 1.0 / safe_y_data

    for i in range(num_points - 1):
        if np.isnan(safe_x_data[i]) or np.isnan(eta_values[i]) or np.isnan(eta_values[i + 1]) \
                or np.isnan(xi_values[i]) or np.isnan(xi_values[i + 1]):
            y_rectifying_vars[i] = np.nan
            continue

        delta_xi = xi_values[i + 1] - xi_values[i]
        if abs(delta_xi) < 1e-12:
            y_rectifying_vars[i] = np.nan
            continue

        k_i = (eta_values[i + 1] - eta_values[i]) / delta_xi
        y_rectifying_vars[i] = k_i * (y_data[i] / x_data[i]) ** 2

    # 5. Вторая разностная производная (y_second_deriv (5))
    # Определена для i = 1, ..., num_points-2
    for i in range(1, num_points - 1):
        y_second_deriv[i] = (y_data[i + 1] - 2 * y_data[i] + y_data[i - 1]) / (h ** 2)

    # Печать результатов
    print(f"{'x':>8} | {'y':>8} | {'y\'_вперед(1)':>15} | {'y\'_центр(2)':>15} | {'y\'_Рунге(3)':>15} | {'y\'_выпр(4)':>15} | {'y\'\'_втор(5)':>15}")
    print("-" * 105)
    for i in range(num_points):
        y_fwd_str = f"{y_fwd_diff[i]:>15.4f}" if not np.isnan(y_fwd_diff[i]) else f"{'Н/Д':>15}"
        y_cen_str = f"{y_central_diff[i]:>15.4f}" if not np.isnan(y_central_diff[i]) else f"{'Н/Д':>15}"
        y_run_str = f"{y_runge[i]:>15.4f}" if not np.isnan(y_runge[i]) else f"{'Н/Д':>15}"
        y_rec_str = f"{y_rectifying_vars[i]:>15.4f}" if not np.isnan(y_rectifying_vars[i]) else f"{'Н/Д':>15}"
        y_sec_str = f"{y_second_deriv[i]:>15.4f}" if not np.isnan(y_second_deriv[i]) else f"{'Н/Д':>15}"

        print(
            f"{x_data[i]:>8.3f} | {y_data[i]:>8.3f} | {y_fwd_str} | {y_cen_str} | {y_run_str} | {y_rec_str} | {y_sec_str}")


# --- Задание 2: Решение ОДУ u'' - 2x^2 u' + 4u = 2x + e^(-x) ---
def task2_solve_ode(n_intervals, alpha_bc, beta_bc, gamma_bc):
    """
    Решает ОДУ методом конечных разностей.
    u'' - 2x^2 u' + 4u = 2x + e^(-x)
    u'(0) = alpha_bc
    u'(1) = beta_bc*u(1) + gamma_bc
    """
    if n_intervals < 2:
        raise ValueError("N_intervals (количество интервалов) должно быть как минимум 2 для ГУ O(h^2).")

    h = 1.0 / n_intervals
    x_grid = np.linspace(0, 1, n_intervals + 1)
    m_points = n_intervals + 1  # Количество точек (уравнений)

    # Коэффициенты для трехдиагональной системы
    a_sys = np.zeros(m_points)  # Поддиагональ
    b_sys = np.zeros(m_points)  # Главная диагональ
    c_sys = np.zeros(m_points)  # Наддиагональ
    d_rhs_sys = np.zeros(m_points)  # Правая часть

    # --- Заполнение коэффициентов для внутренних точек: j = 1, ..., N_intervals-1 ---
    for j in range(1, m_points - 1):
        xj = x_grid[j]
        a_sys[j] = 1 + xj ** 2 * h
        b_sys[j] = -2 + 4 * h ** 2
        c_sys[j] = 1 - xj ** 2 * h
        d_rhs_sys[j] = h ** 2 * (2 * xj + np.exp(-xj))

    # --- Граничное условие при x=0 (j=0) ---
    x1 = x_grid[1]
    coeff_A1_internal = 1 + x1 ** 2 * h
    coeff_B1_internal = -2 + 4 * h ** 2
    coeff_C1_internal = 1 - x1 ** 2 * h
    rhs_D1_internal = h ** 2 * (2 * x1 + np.exp(-x1))

    if abs(coeff_C1_internal) < 1e-14:
        raise ValueError(f"coeff_C1_internal равен нулю (значение: {coeff_C1_internal}) при модификации ГУ при x=0. Попробуйте другое N_intervals.")

    b_sys[0] = coeff_A1_internal - 3 * coeff_C1_internal
    c_sys[0] = coeff_B1_internal + 4 * coeff_C1_internal
    d_rhs_sys[0] = rhs_D1_internal + 2 * h * alpha_bc * coeff_C1_internal

    # --- Граничное условие при x=1 (j=N_intervals) ---
    x_N_minus_1 = x_grid[n_intervals - 1]
    coeff_A_prev_N_internal = 1 + x_N_minus_1 ** 2 * h
    coeff_B_prev_N_internal = -2 + 4 * h ** 2
    coeff_C_prev_N_internal = 1 - x_N_minus_1 ** 2 * h
    rhs_D_prev_N_internal = h ** 2 * (2 * x_N_minus_1 + np.exp(-x_N_minus_1))

    if abs(coeff_A_prev_N_internal) < 1e-14:
        raise ValueError(f"coeff_A_prev_N_internal равен нулю (значение: {coeff_A_prev_N_internal}) при модификации ГУ при x=1. Попробуйте другое N_intervals.")

    N_idx = n_intervals
    a_sys[N_idx] = 4 * coeff_A_prev_N_internal + coeff_B_prev_N_internal
    b_sys[N_idx] = coeff_C_prev_N_internal - (3 - 2 * h * beta_bc) * coeff_A_prev_N_internal
    d_rhs_sys[N_idx] = rhs_D_prev_N_internal - 2 * h * gamma_bc * coeff_A_prev_N_internal

    # Решение трехдиагональной системы
    u_solution = thomas_algorithm(a_sys, b_sys, c_sys, d_rhs_sys)

    return x_grid, u_solution


# --- Основной скрипт выполнения ---
def main():
    print("--- Задание 1: Численное дифференцирование ---")
    # Данные из описания задачи для Задания 1
    x_task1_data = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    y_task1_data = np.array([0.571, 0.889, 1.091, 1.231, 1.333, 1.412])

    try:
        task1_numerical_differentiation(x_task1_data, y_task1_data)
        print("\nАнализ результатов Задания 1 (согласно структуре задачи):")
        print("Столбцы (1)-(4) содержат первые производные, вычисленные разными методами.")
        print("Столбец (5) содержит вторую производную.")
        print(
            "Значения 'Н/Д' (Недоступно/Неприменимо) указывают на точки, где конкретная формула не может быть применена из-за отсутствия соседних точек данных.\n")
    except ValueError as e:
        print(f"Ошибка в Задании 1: {e}")

    print("\n--- Задание 2: Решение ОДУ ---")
    N_ode_intervals = 100  # Количество интервалов (например, 10, 20, 50, 100)
    alpha_ode_bc = 0.1  # Пример значения для u'(0) = alpha
    beta_ode_bc = 0.0  # Пример значения для u'(1) = beta*u(1) + gamma
    gamma_ode_bc = 0.2  # Пример значения для u'(1) = beta*u(1) + gamma

    try:
        x_ode_solution_grid, u_ode_solution = task2_solve_ode(N_ode_intervals, alpha_ode_bc, beta_ode_bc, gamma_ode_bc)

        print(f"ОДУ решено с N={N_ode_intervals} интервалами (h={1.0 / N_ode_intervals:.4f}).")
        print(f"Параметры граничных условий: альфа={alpha_ode_bc}, бета={beta_ode_bc}, гамма={gamma_ode_bc}")
        print("Первые 5 и последние 5 точек решения (x, u(x)):")

        num_solution_points = N_ode_intervals + 1
        indices_to_print = list(range(min(5, num_solution_points))) + \
                           list(range(max(5, num_solution_points - 5), num_solution_points))
        indices_to_print = sorted(list(set(indices_to_print)))

        for i in indices_to_print:
            print(f"x = {x_ode_solution_grid[i]:.3f}, u(x) = {u_ode_solution[i]:.4f}")

        # Построение графика решения
        plt.figure(figsize=(10, 6))
        plt.plot(x_ode_solution_grid, u_ode_solution, label=f'Численное решение (N={N_ode_intervals})')
        plt.title('Решение ОДУ $u\'\' - 2x^2 u\' + 4u = 2x + e^{-x}$')
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.grid(True)
        plt.legend()
        plt.show()

    except ValueError as e:
        print(f"Ошибка в Задании 2: {e}")
    except Exception as e:
        print(f"Произошла непредвиденная ошибка в Задании 2: {e}")

if __name__ == "__main__":
    main()
