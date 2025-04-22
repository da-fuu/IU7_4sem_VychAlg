from printing import *
from task_1d import solve_task_1d, generate_table
from task_2d import solve_task_2d, generate_table_2d
from diff_eq import solve_diff_eq


def main():
    table = generate_table()
    table_2d = generate_table_2d()

    while True:
        print_menu()
        try:
            action = int(input("Выберите действие: "))
        except ValueError:
            print("\nОшибка: ожидался ввод целого числа!")
            continue

        if action == 0:
            break
        elif action == 1:
            print_table(table)
        elif action == 2:
            change_weight(table)
        elif action == 3:
            solve_task_1d(table)
        elif action == 4:
            print_table_2d(table_2d)
        elif action == 5:
            solve_task_2d(table_2d)
        elif action == 6:
            solve_diff_eq()
        elif action > 6:
            print("\nОшибка: ожидался ввод целого числа от 0 до 6!")


if __name__ == "__main__":
    main()
