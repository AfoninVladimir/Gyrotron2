"""v 1.1 введена переменная job

Подпрограмма zbndcrm удаляет столбец с номером col из комплексной матрицы A,
хранимой в ленточном формате в массиве a и возвращает результат
в маccиве b в ленточном формате

Параметры:

  n  - число строк матрицы A (INPUT);
  m  - число столбцов матрицы А (INPUT);
  a - комплексный массив c двойной точностью, в котором хранится разреженная
      матрица A по диагоналям (INPUT)
  lda - размерность массива a по первому индексу, объявленная
        в вызывающей программе (INPUT)
  kla - число ненулевых нижних диагоналей матрицы A (INPUT)
  kua - число ненулевых верхних диагоналей матрицы A (INPUT)
  lowda - номер нижней (левой) диагонали матрицы a (INPUT)
  col - номер удаляемого столбца (INPUT)
job - целая переменная, определяет размещение диагоналей матрицы B в
      массиве b (INPUT):
       job = 0, ненулевые диагонали матрицы сохраняются в массиве b таким образом,
                что правая левая (верхняя диагональ) хранится в строке массива b
                с номером 1, а крайняя левая (нижняя) диагональ хранится в строке
                с номером klb+kub+1. На выходе из программы lowdb =klb+kub+1.
       job = 1, ненулевые диагонали матрицы сохраняются в массиве b таким образом,
                что первые klb строк выделены под ненулевые элементры матрицы b,
                появляющиеся в результате ее разложения методом Гаусса с частичным
                выбором ведущего элемента; на входе в подпрограмму содержимое этих
                строк может быть произвольным. В следующих ku строках с номерами
                от kl+1 до 2*kl+ku+1 хранятся диагонали матрицы b. На выходе из
                программы lowdb =2*klb+kub+1.
       job > 1  диагонали матрицы хранятся в строках  с номерами от lowdw (самая
                левая ненулевая диагональ) до lowdw - (kl+ku+1). Значение lowdw
                при job > 1 задается как входной параметр
  b - комплексный массив c двойной точностью, в котором хранится разреженная
      матрица с удаленным столбцом B по диагоналям (OUTPUT). Число столбцов массива b
      равно m-1
  ldb - размерность массива b по первому индексу, объявленная
        в вызывающей программе (INPUT);
        Размерность массива b по второму индексу должна быть не меньше, чем m
  klb - число ненулевых нижних диагоналей матрицы B (OUT)
  klb - число ненулевых нижних диагоналей матрицы B (OUT)
  kub - число ненулевых верхних диагоналей матрицы B (OUTPUT)
  lowdb - номер нижней (левой) диагонали матрицы b (INPUT)
  ierr - флаг ошибка; если на выходе из подпрограммы ierr.ne.0,
         то произошла ошибка (OUPUT)"""

import numpy as np

zero = 0.0 + 0.0j
one = 1.0 + 0.0j
ImOne = 0.0 + 1.0j


def zbndcrm(n, m, a, lda, kla, kua, lowda, col, job, b, ldb, klb, kub, lowdb, ierr):
    work1 = np.array([[0 + 0j] * len(b[0])] * len(b))

    if col < 1 or col > m:
        ierr = -1
        return

    if lowda > lda or lda < kla + kua + 1 or lowda < kla + kua + 1:
        ierr = -2
        return

    idiaga = lowda - kla

    if col == 1:
        # Удаляется первый столбец
        kub = kua - 1
        klb = kla + 1

        if job == 0:
            lowdb = klb + kub + 1
        if job == 1:
            lowdb = 2 * klb + kub + 1
        idiagb = lowdb - klb
        if lowdb > ldb or ldb < kub + klb + 1 or lowdb < kub + klb + 1:
            ierr = -3
            return

        b[lowdb - (klb + kub + 1):lowdb, 1:m] = zero
        b[lowdb - (klb + kub + 1):lowdb, 1: m - 1] = a[lowda - (kla + kua + 1): lowda, 2: m]
        ierr = 0
        return

    # Удаляется последний столбец

    if col == m:
        kub = kua
        klb = kla
        if job == 0:
            lowdb = klb + kub + 1
        if job == 1:
            lowdb = 2 * klb + kub + 1
        idiagb = lowdb - klb
        b[lowdb - (klb + kub + 1):lowdb, 1:m] = zero
        b[lowdb - (kub + klb + 1): lowdb, 1:m - 1] = a[lowda - (kla + kua + 1): lowda, 1: m - 1]
        ierr = 0
        return

    if 1 < col <= 1 + kua:
        kuw = kua - 1
    if 1 + kua < col < m:
        kuw = kua

    if n >= m + kla:
        klw = kla + 1
    elif n + kua < m and 1 < col < n - kla:
        klw = kla + 1
    elif m > n + kua >= col >= n - kla:
        klw = kla
    elif n + kua < m and n + kua <= col < m:
        klw = kla
    elif n >= m and 1 < col < n - kla:
        klw = kla + 1
    elif n >= m > col >= n - kla:
        klw = kla
    elif n < m and 1 < col < n - kla:
        klw = kla + 1
    elif n < m and n - kla <= col < m:
        klw = kla

    ldw = klw + kuw + 1
    lowdw = ldw
    idiagw = lowdw - klw

    work = np.array([[0 + 0j] * (klw + kuw)] * m)

    work[0: klw + kuw, 0: m-1] = zero
    work = a

    # Копируем верхние диагонали до столбца col-1 включительно
    # for k in range(kua):
    #     for j in range(k, col):
    #         work[idiagw - k, j] = a[idiaga - k, j]

    # Копируем главную и нижние диагонали до столбца col-1 включительно

    # for k in range(kla):
    #     for j in range(col):
    #         work[idiagw + k, j] = a[idiaga + k, j]

    # Копируем верхние диагонали от столбца col+1 до m включительно

    # for k in range(kua):
    #     for j in range(col + 2, min(m, n + k)):
    #         work[idiagw - k + 2, j] = a[idiaga - k, j]

    # Копируем главную и нижнюю диагонали от столбца col+1 до m включительно

    # for k in range(kla):
    #     for j in range(col, min(m, n-k)):
    #         work[idiagw + k + 1, j] = a[idiaga + k, j]

    # Копируем рабочий массив w в выходной массив b

    klb = klw
    kub = kuw
    if job == 0:
        lowdb = klb+kub+1
    if job == 1:
        lowdb = 2 * klb+kub+1
    idiagb = lowdb - klb
    b[0: lowdb, 0: m - 1] = zero

    if ldb < klb + kub + 1 or lowdb > ldb or lowdb < klb+kub+1:
        ierr = -4
        return

    # Верхние диагонали
    for k in range(kuw):
        for j in range(min(m-1, n+k)):
            b[idiagb - k, j] = work[idiagw - k, j]

    # Нижние b и главную
    for k in range(klw):
        for j in range(min(m-1, n-k)):
            b[idiagb + k, j] = work[idiagw + k, j]
