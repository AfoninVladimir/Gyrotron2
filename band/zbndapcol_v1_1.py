"""
Подпрограмма умножает столбец col1 ленточной матрицы aband на скаляр s и
добавляет результат в столбец col2. Результат сохраняется в ленточной
матрице bband

Параметры

n- число строк матриц aband и bband  (INPUT)
m - число столбцов матрицы aband и bband (INPUT)
a - ленточная комплексная матрица с двойной точностью, в которую будет
        добавляться вектор - столбец
lda - длина массива a по первому индексу, объявленная в вызывающей
      программе (INPUT)
kla - число нижних ненулевых диагоналей матрицы a (INPUT)
kua - число верхних ненулевых диагоналей матрицы a (INPUT)
lowda - номер строки с самой нижней (левой) диагональню матрицы a
col1 - номер столбца матрицы a, который будет умножен на скаляр s (INPUT)
col2 - номер столбца матрицы a, в который будет добавлен столбец col1,
       умноженный на скаляр s (INPUT)
s  - действительный скаляр, на который умножается столбец v
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
b - действительный массив для хранения результата сложения в ленточном
    формате (OUT)
ldb - длина массива b по первому индексу, объявленная в вызывающей программе (INPUT)
klb - число нижних ненулевых элементов матрицы b (OUT)
kub - число верхних ненулевых элементов матрицы b (OUT)
lowdb - номер строки с самой нижней (левой) диагональю матрицы b (INPUT)
ierr - флаг ошибки (OUT). Если при выходе из подпрограммы ierr .ne. 0, произошла ошибка"""

zero = 0.0 + 0.0j
one = 1.0 + 0.0j
ImOne = 0.0 + 1.0j
import numpy as np


def zbndapcol(n, m, a, lda, kla, kua, lowda, col1, col2, s, job, b, ldb, klb, kub, lowdb, ierr):
    # Локальные переменные

    work = np.array([0+0j] * (kla + kua + 2))
    work1 = np.array([[0+0j] * len(b[0])] * len(b))

    # Проверка входных параметров
    if kua + kla + 1 > lda:
        ierr = -1
        return ierr

    if lowda > lda:
        ierr = -2
        return ierr

    if (col1 < 1 or col1 > m) or (col2 < 1 or col2 > m):
        ierr = -3
        return ierr

    # Номер строки с главной диагональю матрицы A

    idiaga = lowda - kla

    # Определяем rup и rdown

    rup = max(1, col1 - kua)
    rdown = min(n, col1 + kla)

    # массив для сохранения ненулевых элементов в столбце col1

    work[0:kla + kua + 2] = zero

    for i in range(rup, rdown):
        if work[i - rup + 2] != zero:
            rupnew = i
            break

    # Скорректированные значения rup и rdown

    rupnew = rup
    for i in range(rup, rdown):
        if work[i - rup + 1] != zero:
            rupnew = i
            break

    rdownnew = rdown
    for i in range(rdown, rup, -1):
        if work[i - rup + 1] != zero:
            rdownnew = i
            break

    if rup != rupnew:
        for i in range(rupnew, rdownnew):
            work[i - rupnew + 1] = work[i - rup + 2]

    # значения klb и kub

    kub = kua
    if col2 - rupnew > kua:
        kub = col2 - rupnew

    klb = kla
    if rdownnew - col2 > kla:
        klb = rdownnew - col2

    if job == 0:
        lowdb = klb + kub + 1
    if job == 1:
        lowdb = 2 * klb + kub + 1

    if lowdb > ldb:
        ierr = -4
        return

    if job > 1 and klb + kub + 1 > lowdb:
        ierr = -5
        return

    idiagb = lowdb - klb

    # Копируем массив a в массив b с учетом возможного изменения числа ненулевых диагоналей

    b[idiagb - kub + 1: idiagb + klb + 1, 0:m + 1] = zero
    b[idiagb - kua + 1: idiagb + kla + 1, 0:m + 1] = a[idiaga - kua + 1: idiaga + kla + 1, 0: m + 1]

    # Добавляем вектор work в позиции столбца с номером col

    for i in range(rupnew, rdownnew):
        b[idiagb - col2 + i + 1, col2 + 1] = b[idiagb - col2 + i + 1, col2 + 1] + s * work[i - rupnew + 2]

    # Удаление нулевых внешних диагоналей

    # Нижние дигонали

    kld = 0
    break_for = False
    for k in range(klb + 1, 0, -1):
        if break_for:
            break
        for j in range(min(m, n - k)):
            t = b[idiagb + k + 1, j + 1]
        if t.real != 0 or t.imag != 0:
            klb = k
            break_for = True
            break

    # Верхние диагонали

    kld = 0
    break_for = False
    for k in range(kub, -1, -1):
        for j in range(k + 1, min(m, n + k)):
            t = b[idiagb - k, j]
            if t.real != 0.0 or t.imag != 0.0:
                klb = k
                break_for = True
                break

    # Если число нижних диагоналей скорректировано, опускаем всю ленту
    # диагоналей вниз, чтобы левая (нижняя) диагональ находилась в строке lowdb

    if job == 0:
        lowdbnew = klb + kub + 1
    if job == 1:
        lowdbnew = 2 * klb + kub + 1

    if lowdbnew != lowdb:
        idiagbnew = lowdbnew - klb
        work1[0: klb + kub, 0: m] = b[idiagb - kub: idiagb + klb, 0: m]
        b[idiagbnew - kub - 1: idiagbnew + klb - 1, 0: m] = work1[0: klb + kub, 0: m]

    lowdb = lowdbnew
    ierr = 0

    return
