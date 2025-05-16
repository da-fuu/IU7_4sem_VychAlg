def print_menu():
    print("0. Выйти")
    print("1. Распечатать таблицу")
    print("2. Изменить вес точки")
    print("3. Вывести результаты одномерной аппроксимации")
    print("4. Распечатать таблицу (двумерная аппроксимация)")
    print("5. Вывести результаты двумерной аппроксимации")
    print("6. Решить дифференциальное уравнение")


def print_table(table):
    print("\nСгенерированная таблица\n")
    print("  №    |     X     |     Y    |    W    ")

    for i in range(len(table)):
        print("  {:-3d}  |   {:-5.2f}   |   {:-4.2f}   |   {:-5.2f}   "
              .format(i + 1, table[i][0], table[i][1], table[i][2]))


def print_table_2d(table):
    print("\nСгенерированная таблица (двумерная аппроксимация)\n")
    print("  №    |     X     |     Y     |     Z     |    W    ")

    for i in range(len(table)):
        print("  {:-3d}  |   {:-5.2f}   |   {:-5.2f}   |   {:-5.2f}   |   {:-5.2f}   "
              .format(i + 1, table[i][0], table[i][1], table[i][2], table[i][3]))


def change_weight(table):
    try:
        index = int(input("\nВведите номер точки в таблице: "))
    except ValueError:
        print("\nОшибка: некорректно введён номер точки!")
        return

    if (index > len(table)) or (index < 1):
        print("\nОшибка: в таблице нет точки с таким номером!")
        return

    try:
        weight = float(input("\nВведите новый вес точки: "))
    except ValueError:
        print("\nОшибка: некорректно введён вес точки!")
        return

    table[index - 1][2] = weight
