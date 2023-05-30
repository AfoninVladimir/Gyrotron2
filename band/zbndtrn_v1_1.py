def zbndtrn(n, m, a, lda, maxma, kla, kua, lowda, job, b, ldb, maxmb, klb, kub, lowdb, ierr):
    # Локальные переменные
    zero = 0.0 + 0.0j
    idiaga = 0
    idiagb = 0
    j = 0
    k = 0

    # верхние и нижние диагонали матриц a и b меняются местами

    idiaga = lowda - kla
    klb = kua
    kub = kla

    if job == 0:
        lowdb = klb+kub+1
    if job == 1:
        lowdb = 2 * klb+kub+1
    idiagb = lowdb - klb

    if lda < kla + kua + 1:
        ierr = -1
        return

    if n < 1 or m < 1 or m > maxma:
        ierr = -2
        return

    if lowda > lda or lowdb > ldb:
        ierr = -3
        return

    if lowda < kla+kua+1 or lowdb < klb+kub+1:
        ierr = -4
        return

    # Очищаем массив b (n- число столбцов транспонированной матрицы)

    b[0: ldb, 0: maxmb] = zero

    # Копируем главную и нижние диагонали матрицы A в верхние диагонали матрицы B

    for k in range(kla):
        for j in range(min(m, n-k)):
            b[idiagb - k - 1, j + k] = a[idiaga + k - 1, j]

    # Копируем верхние диагонали матрицы A в нижние диагонали матрицы B

    for k in range(kua):
        for j in range(k+1, min(m, n+k)):
            b[idiagb + k - 1, j - k] = a[idiaga - k - 1, j]

    ierr = 0

    return