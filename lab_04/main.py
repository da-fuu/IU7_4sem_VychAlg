from laplace import find_laplace_argument
from system_solve import solve_system
from boundary import solve_boundary_problem


def main():
    while True:
        print("\n" + "-"*40)
        print("1 - Решить систему уравнений")
        print("2 - Найти аргумент функции Лапласа")
        print("3 - Решить краевую задачу")
        print("0 - Выход")
        print("-"*40)
        
        choice = input("Введите пункт меню: ")
        
        if choice == '1':
            solve_system()
        elif choice == '2':
            find_laplace_argument()
        elif choice == '3':
            solve_boundary_problem()
        elif choice == '0':
            break
        else:
            print("Неверный выбор. Пожалуйста, введите 0, 1, 2 или 3.")


if __name__ == "__main__":
    main()
