import numpy as np
from numpy import *
from scipy.sparse import *


# заменяет "e" на "E" (для Fortran)
def scientific_notation(number, precision_=13):
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


def max_min_comlex(matrix):
    list = [[], []]
    for i in range(len(matrix)):
        max_c = matrix[i][0]
        min_c = matrix[i][0]
        index_max = 0
        index_min = 0
        for j in range(len(matrix[i])):
            if matrix[i][j].real == matrix[i][j].imag == 0:
                continue
            if matrix[i][j].real > max_c.real and matrix[i][j].imag > max_c.imag:
                max_c = matrix[i][j]
                index_max = j
            if matrix[i][j].real < min_c.real and matrix[i][j].imag < min_c.imag:
                min_c = matrix[i][j]
                index_min = j

        list[0].append([f"{scientific_notation(max_c.real, 2)}, {scientific_notation(max_c.imag, 2)}", index_max])
        list[1].append([f"{scientific_notation(min_c.real, 2)}, {scientific_notation(min_c.imag, 2)}", index_min])

    return list


def diff(a, b, src):
    if len(a) != len(b):
        print("Матрицы не совпадают по длине")
        return

    c = []
    c2 = ""
    for i in range(len(a)):
        c.append([])

    accuracy = 1.0e-5
    number_decimal = 2  # кол-во числе после запятой
    for i in range(len(a)):
        if len(a[i]) != len(b[i]):
            print(f"Строки не совпадают по длине: a = {len(a[i])} b = {len(b[i])}")
            return

        for j in range(len(a[i])):
            real_diff = abs(a[i][j].real - b[i][j].real)
            imag_diff = abs(a[i][j].imag - b[i][j].imag)
            c[i].append(complex(real_diff, imag_diff))
            if real_diff > accuracy or imag_diff > accuracy:
                c2 += f"row %2d: col %3d: ({scientific_notation(real_diff, number_decimal)}, {scientific_notation(imag_diff, number_decimal)})\n" % (
                i + 1, j + 1)
    max_min = max_min_comlex(c)

    count = 0
    with open("difference_values_" + src, "w", encoding = "utf-8") as file:

        n = len(str(len(c[0])))
        for i in range(len(c[0])):
            file.write(" " * 8)
            file.write(" " * (14 + n))
            file.write(f"col %{str(n)}d" % (i + 1))
            file.write(" " * (14 - n))
        file.write("\n")

        for i in c:
            count += 1
            file.write(f"row %2d: " % count)
            for j in i:
                file.write(f"({scientific_notation(j.real)}, {scientific_notation(j.imag)}) ")
            file.write("\n")

        file.write("\nMax:\n")
        for j in range(len(max_min[0])):
            if max_min[0][j][0] == "0.00E+00, 0.00E+00":
                continue
            file.write(f"row %2d: col %3d: ({max_min[0][j][0]})\n" % (j + 1, max_min[0][j][1] + 1))

        file.write("\nMin:\n")
        for j in range(len(max_min[1])):
            if max_min[1][j][0] == "0.00E+00, 0.00E+00":
                continue
            file.write(f"row %2d: col %3d: ({max_min[1][j][0]})\n" % (j + 1, max_min[1][j][1] + 1))

        file.write("\nОтносительная разница:\n")
        file.write(c2)

        file.write("\nАбсольтная разница:\n")
        for i in range(len(c)):
            for j in range(len(c[i])):
                if c[i][j] != 0:
                    file.write(
                        f"row %2d: col %3d: ({scientific_notation(c[i][j].real)}, {scientific_notation(c[i][j].imag)}) \n" % (
                            i + 1, j + 1))


def diff_one(a, b, src):
    if len(a) != len(b):
        print("Матрицы не совпадают по длине")
        return

    accuracy = 1.0e-20
    number_decimal = 2  # кол-во числе после запятой
    c = []
    for i in range(len(a)):
        real_diff = abs(a[i].real - b[i].real)
        imag_diff = abs(a[i].imag - b[i].imag)
        c.append(complex(real_diff, imag_diff))
        if real_diff > accuracy or imag_diff > accuracy:
            print(
                f"index %3d: ({scientific_notation(real_diff, number_decimal)}, {scientific_notation(imag_diff, number_decimal)})" % (
                    i))

    with open("difference_values_" + src, "w", encoding = "utf-8") as file:

        n = len(str(len(c)))
        for i in range(len(c)):
            file.write(" " * (17 + n))
            file.write(f"index %{str(n)}d" % (i + 1))
            file.write(" " * (16 - n))
            file.write("|")
        file.write("\n")

        for i in c:
            file.write(f"({scientific_notation(i.real)}, {scientific_notation(i.imag)}) ")
        file.write("\n")
        file.write("\n")

        for i in range(len(c)):
            if c[i] != 0:
                file.write(
                    f"index %3d: ({scientific_notation(c[i].real)}, {scientific_notation(c[i].imag)}) \n" % (i + 1))


########################################################################################################################

# Steps
def steps_diff():
    for i in range(6):
        matrix_1 = []
        matrix_2 = []
        src_matrix = f"MatrixA_v12_Step_{i + 1}.dat"
        print(src_matrix)

        with open(f"{src_matrix}", "r", encoding = "utf-8") as file:
            matrixA = file.readlines()
            for i in range(0, len(matrixA)):
                matrix_1.append([])
                matrixA[i] = matrixA[i].replace("\n", "")
                # print(matrixA[i])
                for j in range(0, len(matrixA[i]), 44):
                    number = matrixA[i][j:j + 43].replace("(", "").replace(")", "").split(",")
                    # print(number)
                    matrix_1[i].append(complex(float(number[0]), float(number[0])))

        with open(f"data/{src_matrix}", "r", encoding = "utf-8") as file:
            matrixA = file.readlines()
            for i in range(0, len(matrixA)):
                matrix_2.append([])
                matrixA[i] = matrixA[i].replace("\n", "")
                for j in range(0, len(matrixA[i]), 44):
                    number = matrixA[i][j:j + 43].replace("(", "").replace(")", "").split(",")
                    matrix_2[i].append(complex(float(number[0]), float(number[0])))

        # Python, Fortran
        diff(matrix_1, matrix_2, src_matrix)


# steps_diff()


# Deltai
def delta_I():
    with open("DELTA_I.dat", "r", encoding = "utf-8") as filePython:
        deltai_py = "".join(filePython.readlines()).replace("(", "").replace(")", "").split("  ")
        for i in range(len(deltai_py)):
            deltai_py[i] = deltai_py[i].split(",")

        deltai_py = deltai_py[0:-2]

    with open("data/DELTA_I.dat", "r", encoding = "utf-8") as fileFortran:
        deltai_for = "".join(fileFortran.readlines()).replace("(", "").replace(")", "").split("  ")
        for i in range(len(deltai_for)):
            deltai_for[i] = deltai_for[i].split(",")
        deltai_for = deltai_for[0:-1]

    for i in range(len(deltai_py)):
        deltai_py[i] = complex(float(deltai_py[i][0]), float(deltai_py[i][1]))
        deltai_for[i] = complex(float(deltai_for[i][0]), float(deltai_for[i][1]))

    diff_one(deltai_py, deltai_for, "DELTA_I.dat")


def Z_I():
    with open("Z_I.dat", "r", encoding = "utf-8") as filePython:
        Z_I_py = "".join(filePython.readlines()).replace("(", "").replace(")", "").split(" ")
    with open("data/Z_I.dat", "r", encoding = "utf-8") as fileFortran:
        Z_I_for = "".join(fileFortran.readlines()).replace("(", "").replace(")", "").split(" ")

    Z_I_py = Z_I_py[0:-1]
    Z_I_for = Z_I_for[0:-1]

    for i in range(len(Z_I_py)):
        Z_I_py[i] = float(Z_I_py[i])
        Z_I_for[i] = float(Z_I_for[i])

    for i in range(len(Z_I_py)):
        if Z_I_py[i] != Z_I_for[i]:
            print(f"index %3d: %.4f => %.4f" % (i, Z_I_py[i], Z_I_for[i]))


def R_I():
    with open("R_I.dat", "r", encoding = "utf-8") as filePython:
        R_I_py = "".join(filePython.readlines()).replace("(", "").replace(")", "").split(" ")
    with open("data/R_I.dat", "r", encoding = "utf-8") as fileFortran:
        R_I_for = "".join(fileFortran.readlines()).replace("(", "").replace(")", "").split(" ")

    R_I_py = R_I_py[0:-1]
    R_I_for = R_I_for[0:-1]

    for i in range(len(R_I_py)):
        R_I_py[i] = float(R_I_py[i])
        R_I_for[i] = float(R_I_for[i])

    accuracy = 1.0e-20
    for i in range(len(R_I_py)):
        dif = abs(R_I_py[i] - R_I_for[i])
        if dif > accuracy:
            print(f"index %3d: {scientific_notation(dif, 2)}" % i)

Z_I()
R_I()