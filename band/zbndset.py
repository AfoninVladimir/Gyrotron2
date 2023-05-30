"""Подпрограмма zbndset присваивает значение val элементу
! комплексной матрицы a с индексами (i,j). Массив a имеет
! двойную точность и хранится в ленточном формате, принятом
! в библиотеке LAPACK
!
! Входные данные
!
!   n - число строк матрицы A
!   m - число столбцов матрицы A
!   a - массив a(lda,*) комплексных чисел с двойной точность
!   lda - размерность массива a по первому индексу, объявленная
!         в вызывающей программе
!   kl - число ненулевых нижних диагоналей матрицы A
!   ku - число ненулевых верхних диагоналей матрицы A
!   lowd - номер строки матрицы a, хранящей нижнюю (левую) диагональ
!   i  - номер строки присваиваемого элемента матрицы
!   j  - номер столбца присваиваемого элемента матрицы
!   val - значение присваиваемого элемента
!
!   ierr - флаг ошибки; если после выхода из программы ierr = 0,
!          ошибки нет
!
!
!----------------------------------------------------------------------*
! Additional details on banded format.  (this closely follows the      *
! format used in linpack. may be useful for converting a matrix into   *
! this storage format in order to use the linpack  banded solvers).    *
!----------------------------------------------------------------------*
!             ---  band storage format  for matrix abd ---             *
! uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of           *
! a in rows of abd starting from the lowest (sub)-diagonal  which  is  *
! stored in row number lowd of abd. the minimum number of rows needed  *
! in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the  *
! j-th  column  of  abd contains the elements of the j-th column of a, *
! from bottom to top: the element a(j+ml,j) is stored in  position     *
! abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.   *
! Generally, the element a(j+k,j) of original matrix a is stored in    *
! position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.               *
! The first dimension nabd of abd must be .ge. lowd                    *
!                                                                      *
!     example [from linpack ]:   if the original matrix is             *
!                                                                      *
!              11 12 13  0  0  0                                       *
!              21 22 23 24  0  0                                       *
!               0 32 33 34 35  0     original banded matrix            *
!               0  0 43 44 45 46                                       *
!               0  0  0 54 55 56                                       *
!               0  0  0  0 65 66                                       *
!                                                                      *
! then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and   *
! if lowd = 5 for example, abd  should be:                             *
!                                                                      *
! untouched --> x  x  x  x  x  x                                       *
!               *  * 13 24 35 46                                       *
!               * 12 23 34 45 56    resulting abd matrix in banded     *
!              11 22 33 44 55 66    format                             *
!  row lowd--> 21 32 43 54 65  *                                       *
!                                                                      *
! * = not used"""


def zbndset(n, m, a, lda, kl, ku, lowd, i, j, val, ierr):
    # Локальные параметры
    # kd - номер строки для храниния главной диагонали
    idiag = 0
    # print(i, j)
    if lda < kl+ku+1:
        # Размерность массива a по первому индексу недостаточна
        # для хранения матрицы с заданным числом верхних и нижних
        # диагоналей
        ierr = -1001
        return ierr

    if i < 0 or i > n or j < 0 or j > m:
        # ошибка в задании индексов вводимого элемента матрицы
        ierr = -1002
        return ierr

    if i > j and j-i > kl:
        # индексы i>j, но элемент лежит вне нижней диагонали
        ierr = -1003
        return ierr

    if j > i and j-i > ku:
        # индексы j>i, но элемент лежит вне верхней диагонали
        ierr = -1004
        return ierr

    if lowd > lda:
        # Номер нижней (левой) диагонали больше размерности
        # массива a по первому индексу
        ierr = -1005
        return ierr

    if lowd < kl+ku+1:
        # Номер нижней (левой) диагонали меньше числа объявленных диагоналей
        # матрицы a
        ierr = -1006
        return ierr

    idiag = lowd - kl

    if i == j:
        # Диагональный элемент
        a[idiag-1, i-1] = val

    elif i > j:
        # элемент на нижней диагонали
        a[idiag + i - j - 1, j-1] = val
    else:
        # на верхней диагонали
        a[idiag + i - j - 1, j-1] = val
    ierr = 0
    return ierr

