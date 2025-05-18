import prettytable as pt


def one_side_left_dif_derivative(h, ys_i, ys_i_l_1=None):
    if ys_i_l_1 is not None:
        return '{0:.3f}'.format((ys_i - ys_i_l_1) / h)
    else:
        return 'None'


def center_dif_derivative(h, ys_p_1=None, ys_m_1=None):
    if ys_p_1 is not None and ys_m_1 is not None:
        return '{0:.3f}'.format((ys_p_1 - ys_m_1) / 2 / h)
    else:
        return 'None'


def second_runge_formula(left_dif_derivative_s_1, left_dif_derivative_s_2, m, p):
    return '{0:.3f}'.format(
        left_dif_derivative_s_1 + (left_dif_derivative_s_1 - left_dif_derivative_s_2) / (m ** p - 1))


def derivative_with_align_vars(theta1, theta2, xsi1, xsi2, y1, x1):
    return '{0:.3f}'.format(((theta2 - theta1) / (xsi2 - xsi1) - 1 / y1) / (- x1 / (y1 ** 2)))


def second_dif_derivative(h, ys_p_1=None, ys=None, ys_m_1=None):
    if ys_p_1 is not None and ys_m_1 is not None:
        return '{0:.3f}'.format((ys_p_1 - 2 * ys + ys_m_1) / (h ** 2))
    else:
        return 'None'


def task1():
    step_1 = 1
    step_2 = 2
    p = 2
    ys = [0.571, 0.889, 1.091, 1.231, 1.333, 1.412]
    xs = [i + 1 for i in range(6)]
    table = pt.PrettyTable()
    columns = ['', 'y', '1', '', '3', '4', '5']

    table.add_column(columns[0], xs)
    table.add_column(columns[1], ys)

    buffer_c = []
    buffer_l = []
    buffer_l_2 = []

    # Разностные
    for i in range(1, len(xs)):
        buffer_l.append(one_side_left_dif_derivative(step_1, ys[i], ys[i - 1]))

    buffer_l.insert(0, one_side_left_dif_derivative(step_1, ys[0]))

    for i in range(1, len(xs) - 1):
        buffer_c.append(center_dif_derivative(step_1, ys[i + 1], ys[i - 1]))

    buffer_c.insert(0, center_dif_derivative(0))
    buffer_c.append('None')

    for i in range(step_2, len(xs)):
        buffer_l_2.append(one_side_left_dif_derivative(step_2, ys[i], ys[i - step_2]))

    # Рунге
    buffer_runge = []

    for i in range(2, len(buffer_l)):
        buffer_runge.append(
            second_runge_formula(float(buffer_l[i]), float(buffer_l_2[i - step_2]), int(step_2 / step_1), p))

    buffer_runge.insert(0, 'None')
    buffer_runge.insert(0, 'None')

    # Выравнивающие
    xsi = [x for x in xs]
    theta = [xs[i] / ys[i] for i in range(len(xs))]

    buffer_align = []

    for i in range(len(xs) - 1):
        buffer_align.append(derivative_with_align_vars(theta[i], theta[i + 1], xsi[i], xsi[i + 1], ys[i], xs[i]))

    buffer_align.append('None')

    # Вторая разностная
    buffer_second = []

    for i in range(1, len(xs) - 1):
        buffer_second.append(second_dif_derivative(step_1, ys[i + 1], ys[i], ys[i - 1]))

    buffer_second.insert(0, center_dif_derivative(0))
    buffer_second.append('None')

    table.add_column('1', buffer_l)
    table.add_column('2', buffer_c)
    table.add_column('3', buffer_runge)
    table.add_column('4', buffer_align)
    table.add_column('5', buffer_second)

    print(table)
