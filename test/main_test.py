from sub_test import sub_test


def fun():
    a = [0, 0, 0]
    ierr = 1
    sub_test(a, ierr)
    print(a, ierr)


fun()

