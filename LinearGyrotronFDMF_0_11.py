import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import sys
import time
from scipy.special import jn as Bessel_JN
from scipy.special import jv as Bessel_JV
from scipy.special import factorial
from scipy.sparse import dia_matrix

from band.zbndset import zbndset
from band.zbndcpy import zbndcpy
from band.zbndtrn_v1_1 import zbndtrn
from band.zbndcrm_v1_1 import zbndcrm
from band.zbndzdrm_v1_1 import zbndzdrm
from band.zbndapcol_v1_1 import zbndapcol
from scipy import special
import numpy as np


def LinearGyrotronFDMFv0_11(iparam, rparam, zi, Ri, Bi):
    #  формируются сразу матрицы в ленточном формате и процедура
    #  исключения столбцов и строк матриц A и B производится
    #  также в ленточном формате

    # Расчетный модуль программы LinearGyrotronFDM, предназначенной нахождениыя комплексных собственных значений
    # и соответствующих собственных векторов линейной теории гиротрона с нефиксированной структурой поля.
    # Input:
    #     iparam(leniparam) - массив с целыми параметрами:
    #     rparam(lenrparam) - массив с действительными параметрами:
    #     zi(lenzi) - массив координат узлов сетки;
    #     Ri(lenRi) - массив значений радиусов резонатора в узлах сетки;
    #     Bi(lenRi) - массив значений магнитного поля в узлах сетки;
    # Output:
    #     vals(lenvals) - массив  собственных чисел;
    #     fields(lenfields) - массив значений собственных функций для полей;
    #     currents(lencurrents) - массив значений собственных функций для токов;
    #     ierr - флаг ошибки

    # Управляющие параметры для процедуры формирования разреженных матриц
    NOUT = 6

    # NEQN - порядок матриц до учета начальных условий для тока
    # NEQN1 - порядок матриц после учета начальных условий
    # ierr - признак ошибки

    neqn = 0
    neqn1 = 0

    # Сдвиг
    omega0 = 0 + 0j
    Shift0 = 0 + 0j

    # Логические флажки

    PrintInfo = False
    PrintDebug = False
    PrintCalcInfo = False
    PrintArraysDebug = False
    ConditionNumberFound = False
    Eigenfunctions = False

    # a, b, fac  - массивы для хранения разреженных матриц в ленточном формате
    # complex list
    a = []
    b = []
    fac = []
    work = []

    # integer
    kla = 1
    kua = 1
    klb = 1
    kub = 1
    kl = 1
    ku = 1
    info = 1
    i = 1
    j = 1
    lworkl = 1
    nev = 1
    ncv = 1
    nconv = 1
    klw = 1
    kuw = 1
    ldw = 1
    lowd = 1
    lowdb = 1
    lowdw = 1
    job = 1

    # Параметры
    zero = complex(0.0, 0.0)
    one = complex(1.0, 0.0)
    ImOne = complex(0.0, 1.0)

    ierr = 1

    # Константы, переменные и массивы для znbandm
    which = 'LM'  # len = 2
    Criteria = ''  # len = 32

    # Логические
    rvec = True
    select = []

    # Целые
    intparam = []
    iwork = []

    # С двойной точностью (real)
    rwork = []
    rd = []

    # Комплексные с двойной точностью

    workd = []
    workl = []
    workev = []
    d = []
    resid = []
    v = []
    ax = []
    bx = []

    # Массивы для параметров рассинхронизма, расстройки магнитного поля и затухания
    # real
    # deltaHi = []
    # zetai = []
    # complex
    # deltai = []
    # di = []

    # Физические константы и масштабирующие множители
    pi = np.pi
    c0 = 299792458.
    e0 = 1.602176620898e-19
    m0 = 9.10938356e-31
    epsilon0 = 8.854187935757958e-12
    mu0 = 4 * pi * 1.0e-7
    mm = 1.0e-3
    GHz = 1.0e+9

    # Целые
    istrF = 1
    istrJ = 1
    istrG = 1

    nrow = 1
    ncol = 1
    ncurr = 1

    # Комплексные
    val = 0 + 0j

    # iparam(leniparam) - массив с целыми параметрами:
    # rparam(lenrparam) - массив с действительными  параметрами:
    # zi(lenzi) - массив координат узлов сетки;
    # Ri(lenRi) - массив значений радиусов резонатора в узлах сетки;
    # Bi(lenRi) - массив значений магнитного поля в узлах сетки;

    # Распаковка параметров

    # Уровень печати модулей пакета ARPACK
    mcaupd = 0
    mceupd = 0

    m = iparam[0]  # - азимутальный индекс
    n = iparam[1]  # - радиальный индекс
    s = iparam[2]  # - номер пространственной гармоники
    NN = iparam[3]  # - число узлов сетки
    NModes = iparam[4]  # - количество рассчитываемых собственных значений и, соответстсвенно, собственных функций

    Eigenfunctions = False  # - признак расчета собственных функций
    if iparam[5] == 1:
        Eigenfunctions = True
    PrintInfo = False  # - признак печати общей информации
    if iparam[6] == 1:
        PrintInfo = True
    PrintDebug = False  # - признак печати отладочной информации
    if iparam[7] == 1:
        PrintDebug = True

    global istart
    istart = iparam[8]  # - номер узла, соответствующий началу пространства взаимодействия
    global ifin
    ifin = iparam[9]  # - номер узла, соответствующий концу пространства взаимодействия

    PrintArraysDebug = False
    if iparam[10] == 1:
        PrintArraysDebug = True

    PrintCalcInfo = False
    if iparam[11] == 1:
        PrintCalcInfo = True

    ConditionNumberFound = False
    if iparam[12] == 1:
        ConditionNumberFound = True

    match iparam[13]:
        case 1:
            Criteria = 'RealPart'
        case 2:
            Criteria = 'ImaginaryPart'
        case 3:
            Criteria = 'Magnitude'
        case 4:
            Criteria = 'ClosestbyRealPart'
        case 5:
            Criteria = 'ClosestbyMagnitude'
        case 6:
            Criteria = 'ClosestbyMImaginaryPart'

    if PrintCalcInfo:
        print("Вход в модуль Python")

    if PrintInfo:
        print("\n" + "Вход в модуль LinearGyrotronFDMF" "\n")

    # ncurr - число узлов сетки, в которых определены значения тока и производной тока
    ncurr = ifin - istart + 1

    # Действительные параметры

    # Re[omega0],     rparam[[11]] - Действительная часть безразмерной частоты линеаризации ГУ 
    # Im[omega0],     rparam[[12]] - Мнимая часть безразмерной частоты линеаризации ГУ 
    # Re[Shift0],     rparam[[13]] - Действительная часть безразмерного сдвига при решении проблемы собственных значений 
    # Im[Shift0],     rparam[[14]] - Мнимая часть безразмерного сдвига при решении проблемы собственных значений 

    R0 = rparam[0]  # - радиус резонатора в плоскости начала выходного рупора
    BB0 = rparam[1]  # - нормировочное значение магнитного поля
    L = rparam[2]  # - длина пространства взаимодействия
    sigma = rparam[3]  # - проводимость стенок резонатора
    zleft = rparam[4]  # - координата начала расчетной области
    zright = rparam[5]  # - координата конца расчетной области
    zs = rparam[6]  # - координата  начала пространства взаимодействия(возможно, скорректированная)
    zf = rparam[7]  # - координата конца пространства взаимодействия(возможно, скорректированная)
    dz = rparam[8]  # - шаг сетки
    tol = rparam[9]  # - погрешность расчета собственных чисел
    omega0 = complex(rparam[10], rparam[11])  # - частота линеаризации ГУ
    Shift0 = complex(rparam[12], rparam[13])  # - сдвиг при решении обобщенной проблемы собстрвенных значений
    VBeam = rparam[14]  # - Напряжение пучка
    IBeam = rparam[15]  # - Ток пучка
    RBeam = rparam[16]  # - радиус пучка
    g = rparam[17]  # - питч - фактор
    numn = rparam[18]  # - корень производной функции Бесселя
    ConductanceThreshold = rparam[19]  # - порог учета проводимости стенок резонатора

    if PrintCalcInfo:
        print("Массивы iparam и rparam декодированы")

    # Выделяем память под массивы
    di = []  # len = NN
    deltai = []  # len = NN
    DeltaHi = []  # len = NN
    zetai = []  # len = NN

    if PrintCalcInfo:
        print("Массивы di, deltai, deltaHi, zetai выделены")

    # Критические частоты (круговая и циклическая)
    omegamn = c0 * numn / R0
    f0 = omegamn / (2 * pi)

    # Релятивистские факторы
    gamma0 = 1 + e0 * VBeam / (m0 * c0 ** 2)
    beta0 = np.sqrt(1 - 1 / (gamma0 ** 2))
    betalong = beta0 / np.sqrt(1 + g ** 2)
    betatrans = g * betalong

    # Безразмерный параметр численной схемы
    alpha = complex(4 * s / (g ** 2), 0.0)

    # Циклотронные    частоты
    # gyroratio - гиромагнитное    отношени

    gyroratio = e0 / (m0 * gamma0)
    OmegaH = gyroratio * BB0
    fH = OmegaH / (2 * pi)

    # Параметр расстройки критической частоты от частоты циклотронн гармоники
    DeltaH = 2 / (betatrans**2) * (omegamn - s * OmegaH) / OmegaH

    # Нормировочный множитель для координаты
    znorm = 2 * betalong/betatrans**2 * c0/OmegaH

    # Сетка
    zetai = zi / znorm
    dzeta = dz / znorm

    # Параметры    численной    схемы
    dzeta2 = dzeta * dzeta
    alpha = complex(4 * s / (g ** 2), 0.0)
    alphadiv2dzeta = alpha / (2 * dzeta)
    onediv2dzeta = complex(1.0 / (2 * dzeta), 0.0)
    onedivdzeta2 = complex(1.0 / dzeta2, 0.0)
    twodivdzeta = complex(2.0 / dzeta, 0.0)
    twodivdzeta2 = complex(2.0 / dzeta2, 0.0)

    if PrintCalcInfo:
        print("Вычисленная сетка конечных разностей")

    # Безразмерные параметры взаимодействия и тока
    if m == 0:
        ms = abs(s)
    elif m > 0:
        ms = abs(m) - s
    else:
        ms = abs(m) + s

    G0 = Bessel_JN(ms, (numn * RBeam / R0)) ** 2 / ((numn ** 2 - m ** 2) * Bessel_JN(abs(m), numn) ** 2)
    I0 = 64 * (e0 * IBeam) / (m0 * c0) * mu0 / (4 * pi) * betalong * (betatrans ** (2 * (s - 4))) / gamma0 * (
            s ** 4) * (float((s ** s)) / float(2 ** s * factorial(s))) ** 2 * G0

    # Массивы расстроек радиуса резонатора и магнитного поля
    for i in range(len(Ri)):
        deltai.append(complex(((4 * (s ** 2) * (betalong ** 2) / (betatrans ** 4)) * (1.0 - ((R0 / Ri[i]) ** 2))), 0.0))
    DeltaHi = (2 / betatrans ** 2) * (omegamn - s * gyroratio * Bi) / OmegaH

    # Расчет массива параметров затухания
    # deltas - толщина скин - слоя

    if sigma > ConductanceThreshold:
        deltas = 0.0
        for i in range(NN):
            di.append(zero)
    else:
        deltas = np.sqrt(2 / (sigma * mu0 * omegamn))
        for i in range(NN):
            di.append((1.0 + ImOne) * znorm ** 2 * (numn / Ri[i]) ** 2 * (deltas / Ri[i]) *
                      (1.0 + m ** 2 / (numn ** 2 - m ** 2) * (omegamn * Ri[i] / (numn * c0)) ** 2))

    # deltai -= di
    for i in range(len(deltai)):
        deltai[i] -= di[i]

    #############
    # with open("DELTA_I.dat", "w", encoding = "utf-8") as file:
    #     for i in deltai:
    #         if str(i.real)[0] == "-" and str(i.imag)[0] == "-":
    #             file.write(f"({scientific_notation(i.real)},{scientific_notation(i.imag)})")
    #         elif str(i.real)[0] != "-" and str(i.imag)[0] == "-":
    #             file.write(f"( {scientific_notation(i.real)},{scientific_notation(i.imag)})")
    #         elif str(i.real)[0] == "-" and str(i.imag)[0] != "-":
    #             file.write(f"({scientific_notation(i.real)}, {scientific_notation(i.imag)})")
    #         else:
    #             file.write(f"( {scientific_notation(i.real)}, {scientific_notation(i.imag)})")
    #         file.write("  ")
    #############

    if PrintCalcInfo:
        print("Массивы deltai и DeltaHi вычислены")

    # Печать информации в консольное окно фортрана
    if PrintInfo:
        print("Входные параметры:")
        print(f"Азимутальный индекс m = {m}")
        print(f"Радиальный индекс n = {n}")
        print(f"Гармоническое число циклотрона s = {s}")
        print(f"Радиус однородной части резонатора R0 = {R0 / mm:12.5f} мм")
        print(f"Длина однородной части резонатора L =  {L / mm:12.5f} мм")
        print(f"Координата начальной точки области моделирования zleft = {zleft / mm} мм")
        print(f"Координата конечной точки области моделирования zright = {zright / mm} мм")
        print(f"Координата начальной точки области взаимодействия zstart = {zs / mm} мм")
        print(f"Координата конечной точки области взаимодействия zfin = {zf / mm} мм")
        print(f"Шаг сетки dz = {dz / mm} мм")
        print(f"Количество рассчитанных режимов NModes = {NModes}")
        print(f"Напряжение пучка VBeam = {VBeam} V")
        print(f"Ток пучка iBeam = {IBeam} A")
        print(f"Радиус пучка RBeam = {RBeam / mm :12.5f} мм")
        print(f"Коэффициент шага электронного пучка g = {g:12.5f}")
        print(f"Нормировочные константы для магнитного поля B0 = {BB0:12.5f} T")
        if sigma > ConductanceThreshold:
            print("Идеально проводящие стенки резонатора")
        else:
            print(f"Проводимость стенок резонатора sigma = {sigma} См/м")
        print("Рассчитанные параметры:")
        print(f"Корень производной функции Бесселя num = {numn: 14.7f}")
        print(f"Частота среза режима TEmn f = {f0 / GHz : 12.5f} ГГц")
        print(f"Релятивистский фактор gamma = {gamma0}")
        print(f"Релятивистский фактор beta = {beta0}")
        print(f"Продольный релятивистский фактор betalon = {betalong}")
        print(f"Трансверсальный релятивистский фактор betatran = {betatrans}")
        print(f"Циклотронная частота f = {fH / GHz : 12.5f} ГГц")
        print(f"Несоответствие циклотронного резонанса Delta = {DeltaH}")
        print(f"Безразмерный параметр поля пучка G = {G0}")
        print(f"Безразмерный текущий параметр I = {I0}")
        print(f"Константа нормализации для координаты znorm = {znorm / mm : 12.5f}, мм")
        print(f"Количество точек сетки NN = {NN}")
        print(f"Номер первой точки области взаимодействия istart = {istart}")
        print(f"Номер последней точки области взаимодействия ifin   = {ifin}")
        # print(f"Максимальное число ненулевых элементов в разреженной матрице MAXNZ = {MAXNZ}")
        print(f"Безразмерный размер сетки dzeta = {dzeta}")
        print(f"")

    """ Формирование матриц """
    # neqn - размерность матриц до удаления столбцов и строк, проводимого для учета
    # граничных условий для тока и его производной в начале пространства взаимодействия
    neqn = indexF(NN)

    """ Определение числа верхних и нижних диагоналей матрицы A """
    kla = 8
    kua = 3

    lda = 2 * kla + kua + 1
    lowda = kla + kua + 1
    lda = 12

    """Формирование матрицы A"""
    print(lda, neqn)
    a = np.zeros((lda, neqn), dtype = complex)

    # Формирование коэффициентов уравнений для поля:

    # Диагональные элементы:
    val = 0.5 * (twodivdzeta2 - deltai[0]) - 0.5 * omega0 + ImOne * np.sqrt(deltai[0] + omega0) / dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexF(1), indexF(1), val, ierr)
    if ierr != 0:
        print("Error 404:", ierr)
        return

    val = 0.5 * (twodivdzeta2 - deltai[NN - 1]) - 0.5 * omega0 + ImOne * np.sqrt(deltai[NN - 1] + omega0) / dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexF(NN), indexF(NN), val, ierr)
    if ierr != 0:
        print("Error 410:", ierr)
        return

    for i in range(1, NN - 1):
        val = twodivdzeta2 - deltai[i] - omega0
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexF(i + 1), indexF(i + 1), val, ierr)
    if ierr != 0:
        print("Error 417:", ierr)
        return

    with open("MatrixA_v12_Step_1.dat", "w", encoding = "utf-8") as matrixA:
        for i in a:
            for j in i:
                if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                    matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                else:
                    matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                matrixA.write(" ")
            matrixA.write("\n")
    print("MatrixA_v12_Step_1.dat")

    # Элементы нижней диагонали
    val = -onedivdzeta2
    for i in range(1, NN):
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexF(i + 1), indexF(i), val, ierr)
    if ierr != 0:
        print("Error 425:", ierr)
        return

    # Элементы верхней диагонали
    val = -onedivdzeta2
    for i in range(0, NN - 1):
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexF(i + 1), indexF(i + 2), val, ierr)
        if ierr != 0:
            print("Error 436:", ierr)
            return

    with open("MatrixA_v12_Step_2.dat", "w", encoding = "utf-8") as matrixA:
        for i in a:
            for j in i:
                if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                    matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                else:
                    matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                matrixA.write(" ")
            matrixA.write("\n")
    print("MatrixA_v12_Step_2.dat")

    # Элементы, умножаемые на ток
    val = complex(I0, 0.0)
    for i in range(istart, ifin):
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexF(i), indexJ(i), val, ierr)
        if ierr != 0:
            print("Error 443:", ierr)
            return

    with open("MatrixA_v12_Step_3.dat", "w", encoding = "utf-8") as matrixA:
        for i in a:
            for j in i:
                if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                    matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                else:
                    matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                matrixA.write(" ")
            matrixA.write("\n")
    print("MatrixA_v12_Step_3.dat")

    # Формирование коэффициентов уравнения для тока
    for i in range(istart + 1, ifin - 1):

        # Диагональные элементы
        val = ImOne * alpha * DeltaHi[i] + ImOne * omega0
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(i), indexJ(i), val, ierr)
        if ierr != 0:
            print("Error 454:", ierr)
            return

        # Элементы нижней диагонали
        val = -alphadiv2dzeta
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(i), indexJ(i - 1), val, ierr)
        if ierr != 0:
            print("Error 461:", ierr)
            return

        # Элементы верхней диагонали
        val = alphadiv2dzeta
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(i), indexJ(i + 1), val, ierr)
        if ierr != 0:
            print("Error 468:", ierr)
            return

        # Коэффициент, умножаемый на производную тока
        val = -alpha
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(i), indexG(i), val, ierr)
        if ierr != 0:
            print("Error 475:", ierr)
            return

    with open("MatrixA_v12_Step_4.dat", "w", encoding = "utf-8") as matrixA:
        for i in a:
            for j in i:
                if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                    matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                else:
                    matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                matrixA.write(" ")
            matrixA.write("\n")
    print("MatrixA_v12_Step_4.dat")

    # Коэффициенты в строке indexJ(ifin)
    val = 3 * alphadiv2dzeta + ImOne * alpha * DeltaHi[ifin] + ImOne * omega0
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(ifin), indexJ(ifin), val, ierr)
    if ierr != 0:
        print("Error 481:", ierr)
        return

    val = -4 * alphadiv2dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(ifin), indexJ(ifin - 1), val, ierr)
    if ierr != 0:
        print("Error 488:", ierr)
        return

    val = alphadiv2dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(ifin), indexJ(ifin - 2), val, ierr)
    if ierr != 0:
        print("Error 494:", ierr)
        return

    # Коэффициент, умножаемый на производную тока в этой строке
    val = -alpha
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexJ(ifin), indexG(ifin), val, ierr)
    if ierr != 0:
        print("Error 501:", ierr)
        return

    with open("MatrixA_v12_Step_5.dat", "w", encoding = "utf-8") as matrixA:
        for i in a:
            for j in i:
                if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                    matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                else:
                    matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                matrixA.write(" ")
            matrixA.write("\n")
    print("MatrixA_v12_Step_5.dat")

    # Коэффициенты уравнений для производной тока
    for i in range(istart + 1, ifin - 1):

        # Диагональные элементы
        val = ImOne * alpha * DeltaHi[i] + ImOne * omega0
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(i), indexG(i), val, ierr)
        if ierr != 0:
            print("Error 512:", ierr)
            return

        # Элементы ниже главной диагонали
        val = -alphadiv2dzeta
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(i), indexG(i - 1), val, ierr)
        if ierr != 0:
            print("Error 519:", ierr)
            return

        # Элементы выше главной диагонали
        val = alphadiv2dzeta
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(i), indexG(i + 1), val, ierr)
        if ierr != 0:
            print("Error 526:", ierr)
            return

        # Элементы, умножаемые на F(i)
        val = ImOne * alphadiv2dzeta
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(i), indexF(i - 1), val, ierr)
        if ierr != 0:
            print("Error 533:", ierr)
            return

        # Элементы, умножаемые на F(i-1)
        val = ImOne * alphadiv2dzeta
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(i), indexF(i - 1), val, ierr)
        if ierr != 0:
            print("Error 540:", ierr)
            return

        # Элементы, умножаемые на F(i+1)
        val = - ImOne * alphadiv2dzeta
        ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(i), indexF(i + 1), val, ierr)
        if ierr != 0:
            print("Error 547:", ierr)
            return

    # Строка indexG(ifin)
    val = -3 * ImOne * alphadiv2dzeta + alpha * (DeltaHi[ifin] - 1.0) + omega0
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(ifin), indexF(ifin), val, ierr)
    if ierr != 0:
        print("Error 554:", ierr)
        return

    val = 4 * ImOne * alphadiv2dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(ifin), indexF(ifin - 1), val, ierr)
    if ierr != 0:
        print("Error 560:", ierr)
        return

    val = -ImOne * alphadiv2dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(ifin), indexF(ifin - 2), val, ierr)
    if ierr != 0:
        print("Error 566:", ierr)
        return

    val = 3 * alphadiv2dzeta + ImOne * alpha * DeltaHi[ifin]
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(ifin), indexG(ifin), val, ierr)
    if ierr != 0:
        print("Error 572:", ierr)
        return

    val = -4 * alphadiv2dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(ifin), indexG(ifin - 1), val, ierr)
    if ierr != 0:
        print("Error 578:", ierr)
        return

    val = alphadiv2dzeta
    ierr = zbndset(neqn, neqn, a, lda, kla, kua, lowda, indexG(ifin), indexG(ifin - 2), val, ierr)
    if ierr != 0:
        print("Error 584:", ierr)
        return

    if PrintCalcInfo:
        pass
        # print("Ленточная матрица А сформирована")

    # Матрица A сформирована
    with open("MatrixA_v12_Step_6.dat", "w", encoding = "utf-8") as matrixA:
        for i in a:
            for j in i:
                if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                    matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                    matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                else:
                    matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                matrixA.write(" ")
            matrixA.write("\n")
    print("MatrixA_v12_Step_6.dat")
    # print("Матрица A сформирована")
    return

    """Формирование матрицы B"""

    klb = 2
    kub = 0
    ldb = lda
    lowdb = klb + kub + 1

    if PrintDebug:
        print()
        print(f"Количество нижних диагоналей матрицы B klb = {klb}")
        print(f"Количество верхних диагоналей матрицы B kub = {kub}")
        print()

    # Выделение памяти под матрицу b
    b = np.zeros((lda, neqn), dtype = complex)

    # Коэффициенты уравнений для поля
    val = 0.5 * (1.0 - ImOne / (dzeta * np.sqrt(deltai[1] + omega0)))
    ierr = zbndset(neqn, neqn, b, ldb, klb, kub, lowdb, indexF(1), indexF(1), val, ierr)
    if ierr != 0:
        print("Error 600:", ierr)
        return

    val = 0.5 * (1.0 - ImOne / (dzeta * np.sqrt(deltai[NN - 1] + omega0)))
    ierr = zbndset(neqn, neqn, b, ldb, klb, kub, lowdb, indexF(NN), indexF(NN), val, ierr)

    for i in range(2, NN - 1):
        val = 1.0 + 0.0j
        ierr = zbndset(neqn, neqn, b, ldb, klb, kub, lowdb, indexF(i), indexF(i), val, ierr)
        if ierr != 0:
            print("Error 610:", ierr)
            return

    # Коэффициенты в уравнениях для тока и производной тока
    for i in range(istart + 1, ifin + 1):
        # Для тока
        val = -ImOne
        ierr = zbndset(neqn, neqn, b, ldb, klb, kub, lowdb, indexJ(i), indexJ(i), val, ierr)
        if ierr != 0:
            print("Error 619:", ierr)
            return

        # Для производной тока
        val = -ImOne
        ierr = zbndset(neqn, neqn, b, ldb, klb, kub, lowdb, indexG(i), indexG(i), val, ierr)
        if ierr != 0:
            print("Error 626:", ierr)
            return

        val = -1.0 + 0.0j
        ierr = zbndset(neqn, neqn, b, ldb, klb, kub, lowdb, indexG(i), indexF(i), val, ierr)
        if ierr != 0:
            print("Error 632:", ierr)
            return

    if PrintCalcInfo:
        print("Ленточная матрица B сформирована")

    # Матрица B сформирована

    if PrintArraysDebug:
        # Вывод матриц A и B до учета граничных условий
        with open("MatrixA_v12_1.dat", "w", encoding = "utf-8") as matrixA:
            for i in a:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixA.write(" ")
                matrixA.write("\n")

        with open("MatrixB_v12_1.dat", "w", encoding = "utf-8") as matrixB:
            for i in b:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixB.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixB.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixB.write(" ")
                matrixB.write("\n")

    """ Учет граничныx условий для производной тока  """

    istrF = indexF(istart)
    istrG = indexG(istart)
    istrJ = indexJ(istart)

    # work - рабочий массив

    ldw = lda
    work = np.zeros((ldw, neqn), dtype = complex)

    # Извлечение из матрицы A столбца номер istrG, умножение на ImOne и добавление в столбец istrF

    job = 0
    zbndapcol(neqn, neqn, a, lda, kla, kua, lowda, istrG, istrF, ImOne, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print("Error 671")
        print("Ошибка в zbndapcol при преобразовании матрицы A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Сохраняем результат из матрицы work в матрицу A
    zbndcpy(neqn, neqn, work, ldw, klw, kuw, lowdw, job, a, lda, kla, kua, lowda, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при копировании массива work в массив A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Извлечение из матрицы B столбца номер istrG, умножение на ImOne и добавление в столбец istrF
    zbndapcol(neqn, neqn, b, ldb, klb, kub, lowdb, istrG, istrF, ImOne, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print("Ошибка в zbndapcol при преобразовании матрицы B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Сохраняем результат из матрицы work в матрицу B
    zbndcpy(neqn, neqn, work, ldw, klw, kuw, lowdw, job, b, ldb, klb, kub, lowdb, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при копировании массива work в массив B")
        print(f"Код ошибки ierr = {ierr}")
        return

    if PrintArraysDebug:
        # Вывод в файлы матрицы A и B после учета граничных условий для производной тока
        with open("MatrixA_v12_2.dat", "w", encoding = "utf-8") as matrixA:
            for i in a:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixA.write(" ")
                matrixA.write("\n")

        with open("MatrixB_v12_2.dat", "w", encoding = "utf-8") as matrixB:
            for i in b:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixB.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixB.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixB.write(" ")
                matrixB.write("\n")

    """Удаляем строки и столбцы из матриц A и B для учета ГУ для производной поля"""
    nrow = neqn
    ncol = neqn

    # Матрица A

    # Удаляем из матрицы A столбец istrG
    job = 0
    zbndcrm(nrow, ncol, a, lda, kla, kua, lowda, istrG, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrG из матрицы A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем из матрицы A столбец istrJ
    zbndcrm(nrow, ncol - 1, work, ldw, klw, kuw, lowdw, istrJ, job, a, lda, kla, kua, lowda, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrJ из матрицы A")
        print(f"Код ошибки ierr = {ierr}")
        return

    #  Транспонирование матрицы  A
    zbndtrn(nrow, ncol - 2, a, lda, neqn, kla, kua, lowda, job, work, ldw, neqn, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndtrn при первом транспонировании массива A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем из матрицы transp(A) столбец istrG
    zbndcrm(ncol - 2, nrow, work, ldw, klw, kuw, lowdw, istrG, job, a, lda, kla, kua, lowda, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrG из матрицы A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем из матрицы transp(A) столбец istrJ
    zbndcrm(nrow - 2, nrow - 1, a, lda, kla, kua, lowda, istrJ, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrJ из транспонированной матрицы A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Обратное транспонирование матрицы A
    zbndtrn(nrow - 2, ncol - 2, work, ldw, neqn, klw, kuw, lowdw, job, a, lda, neqn, kla, kua, lowda, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndtrn при втором транспонировании массива A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем нулевые нижние и верхние диагонали матрицы A
    zbndzdrm(neqn1, neqn1, a, lda, kla, kua, lowda, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndzdrm при удалении пустых диагоналей из массива A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Копируем в матрицу A
    job = 0
    zbndcpy(neqn1, neqn1, work, ldw, klw, kuw, lowdw, job, a, lda, kla, kua, lowda, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при копировании рабочего массива work в массив A")
        print(f"Код ошибки ierr = {ierr}")
        return

    if PrintDebug:
        print(f"Количество нижних и верхних диагоналей матрицы A"
              f" после удаления пустых диагоналей с учетом равно (kla,kua) = ({kla},{kua})")

    nrow = neqn
    ncol = neqn

    """Матрица B"""

    # Удаляем из матрицы B столбец istrG
    job = 0
    zbndcrm(nrow, ncol, b, ldb, klb, kub, lowdb, istrG, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrG из матрицы B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем из матрицы B столбец istrJ
    zbndcrm(nrow, ncol - 1, work, ldw, klw, kuw, lowdw, istrJ, job, b, ldb, klb, kub, lowdb, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrJ из матрицы B")
        print(f"Код ошибки ierr = {ierr}")
        return

    #  Транспонирование матрицы  B
    zbndtrn(nrow, ncol - 2, b, ldb, neqn, klb, kub, lowdb, job, work, ldw, neqn, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndtrn при первом транспонировании массива B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем из матрицы transp(B) столбец istrG
    zbndcrm(ncol - 2, nrow, work, ldw, klw, kuw, lowdw, istrG, job, b, ldb, klb, kub, lowdb, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrG из матрицы B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем из матрицы transp(B) столбец istrJ
    zbndcrm(nrow - 2, nrow - 1, b, ldb, klb, kub, lowdb, istrJ, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcrm при удалении строки istrJ из транспонированной матрицы B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Обратное транспонирование матрицы b
    zbndtrn(nrow - 2, ncol - 2, work, ldw, neqn, klw, kuw, lowdw, job, b, ldb, neqn, klb, kub, lowdb, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndtrn при втором транспонировании массива B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Удаляем нулевые нижние и верхние диагонали матрицы B
    zbndzdrm(neqn1, neqn1, b, ldb, klb, kub, lowdb, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndzdrm при удалении пустых диагоналей из массива B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Копируем в матрицу A
    job = 0
    zbndcpy(neqn1, neqn1, work, ldw, klw, kuw, lowdw, job, b, ldb, klb, kub, lowdb, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при копировании рабочего массива work в массив B")
        print(f"Код ошибки ierr = {ierr}")
        return

    if PrintDebug:
        print(f"Количество нижних и верхних диагоналей матрицы B"
              f" после удаления пустых диагоналей с учетом равно (klb,kub) = ({klb},{kub})")

    # Вывод матриц в файл после удаления строк и столбцов

    if PrintArraysDebug:
        # Вывод матриц A и B до учета граничных условий
        with open("MatrixA_v12_3.dat", "w", encoding = "utf-8") as matrixA:
            for i in a:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixA.write(" ")
                matrixA.write("\n")

        with open("MatrixB_v12_3.dat", "w", encoding = "utf-8") as matrixB:
            for i in b:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixB.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixB.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixB.write(" ")
                matrixB.write("\n")

    # Сдвигаем ленты диагоналей матриц A и B таким образом, чтбы оснободить место по появляющиеся
    # ненулевые диагонали в результате разложения Гаусса.

    # Массив a
    job = 1
    zbndcpy(neqn1, neqn1, a, lda, kla, kua, lowda, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при смещении диагоналей матрицы A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Полная очистка массива a
    a[0:lda, 0:neqn] = zero
    zbndcpy(neqn1, neqn1, work, ldw, klw, kuw, lowdw, job, a, lda, kla, kua, lowda, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при копировании сдвинутых диагоналей матрицы A")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Массив b
    job = 2
    lowdw = lowda - kla + klb
    work[0:ldw, 0:neqn] = zero
    zbndcpy(neqn1, neqn1, b, ldb, klb, kub, lowdb, job, work, ldw, klw, kuw, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при смещении диагоналей матрицы B")
        print(f"Код ошибки ierr = {ierr}")
        return

    # Полная очистка массива b
    b[0:ldb, 0:neqn] = zero
    job = 2
    lowdb = lowdw
    zbndcpy(neqn1, neqn1, work, ldw, klw, kuw, lowdw, job, b, ldb, klb, kub, lowdw, ierr)
    if ierr != 0:
        print()
        print("Ошибка в zbndcpy при копировании сдвинутых диагоналей матрицы B")
        print(f"Код ошибки ierr = {ierr}")
        return

    klb = kla
    kub = kua
    lowdb = lowda

    if PrintArraysDebug:
        # Вывод матриц A и B в файлы на диске перед вызовом программы решения матричной проблемы собственных значений
        with open("MatrixA_v12_4.dat", "w", encoding = "utf-8") as matrixA:
            for i in a:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixA.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixA.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixA.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixA.write(" ")
                matrixA.write("\n")

        with open("MatrixB_v12_4.dat", "w", encoding = "utf-8") as matrixB:
            for i in b:
                for j in i:
                    if str(j.real)[0] == "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"({scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] != "-" and str(j.imag)[0] == "-":
                        matrixB.write(f"( {scientific_notation(j.real)},{scientific_notation(j.imag)})")
                    elif str(j.real)[0] == "-" and str(j.imag)[0] != "-":
                        matrixB.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    else:
                        matrixB.write(f"( {scientific_notation(j.real)}, {scientific_notation(j.imag)})")
                    matrixB.write(" ")
                matrixB.write("\n")

    #  Матрицы A и B полностью сформированы для решения матричной проблемы собственных значений
    # (Сдвиг на Shift0 проводится в подпрограмме znbandm)

    if PrintCalcInfo:
        print("\n" + "Матрицы A и B сформированы полностью")

    """Решение обобщенной проблемы собственных  значений """


def indexF(i):
    # Номер строки и столбца переменной поля в узле i

    if i <= istart:
        return i
    elif istart < i <= ifin:
        return istart + 3 * (i - istart)
    elif ifin < i:
        return istart + 3 * (ifin + 1 - istart) + (i - ifin - 1)
    else:
        return -1


def indexJ(i):
    # Номер строки и столбца переменной поля в узле i
    if i < istart:
        return 0
    elif istart <= i <= ifin:
        return indexF(i) + 1
    elif ifin < i:
        return 0
    else:
        return -1


def indexG(i):
    # Номер строки и столбца переменной поля в узле i
    if i < istart:
        return 0
    elif istart <= i <= ifin:
        return indexF(i) + 2
    elif ifin < i:
        return 0
    else:
        return -1


# заменяет "e" на "E" (для Fortran)
def scientific_notation(number, precision_ = 13):
    """
    Format a floating-point scalar as a decimal string in scientific notation.

    Provides control over rounding, trimming and padding. Uses and assumes
    IEEE unbiased rounding. Uses the "Dragon4" algorithm.

    Parameters
    ----------
    x : python float or numpy floating scalar
        Value to format.
    precision : non-negative integer or None, optional
        Maximum number of digits to print. May be None if `unique` is
        `True`, but must be an integer if unique is `False`.
    unique : boolean, optional
        If `True`, use a digit-generation strategy which gives the shortest
        representation which uniquely identifies the floating-point number from
        other values of the same type, by judicious rounding. If `precision`
        was omitted, print all necessary digits, otherwise digit generation is
        cut off after `precision` digits and the remaining value is rounded.
        If `False`, digits are generated as if printing an infinite-precision
        value and stopping after `precision` digits, rounding the remaining
        value.
    trim : one of 'k', '.', '0', '-', optional
        Controls post-processing trimming of trailing digits, as follows:

        * 'k' : keep trailing zeros, keep decimal point (no trimming)
        * '.' : trim all trailing zeros, leave decimal point
        * '0' : trim all but the zero before the decimal point. Insert the
          zero if it is missing.
        * '-' : trim trailing zeros and any trailing decimal point
    sign : boolean, optional
        Whether to show the sign for positive values.
    pad_left : non-negative integer, optional
        Pad the left side of the string with whitespace until at least that
        many characters are to the left of the decimal point.
    exp_digits : non-negative integer, optional
        Pad the exponent with zeros until it contains at least this many digits.
        If omitted, the exponent will be at least 2 digits."""
    return str(np.format_float_scientific(number, unique = False, precision = precision_, trim = 'k',
                                          exp_digits = 2)).replace("e", "E")
