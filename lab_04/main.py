from laplace import find_laplace_argument
from system_solve import solve_system
from boundary import solve_boundary_problem


def main():
    while True:
        print("\n" + "="*40)
        print(f"{'МЕНЮ':^40}")
        print("="*40)
        print("1 - Решить систему уравнений")
        print("2 - Найти аргумент функции Лапласа")
        print("3 - Решить краевую задачу")
        print("0 - Выход")
        print("="*40)
        
        choice = input("Введите пункт меню: ")
        
        if choice == '1':
            print("\nРешение системы уравнений...")
            solve_system()
        elif choice == '2':
            print("\nПоиск аргументов функции Лапласа...")
            find_laplace_argument()
        elif choice == '3':
            print("\nРешение краевой задачи...")
            solve_boundary_problem()
        elif choice == '0':
            print("Выход...")
            break
        else:
            print("Неверный выбор. Пожалуйста, введите 0, 1, 2 или 3.")


if __name__ == "__main__":
    main()
