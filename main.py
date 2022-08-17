import numpy as np


def generate_polynom(number, e):
    number_bin = bin(number)[2:]
    # return np.array([int(d) for d in number_bin.zfill(e)])
    return np.array([int(d) for d in number_bin])


def generate_sx_polynom(diff):
    s = "1"
    for i in range(diff):
        s += "0"
    return np.array([int(d) for d in s])


def calculate_mod(polynom_multiplication, irreducible_polynom, e):
    while True:
        if len(polynom_multiplication) < len(irreducible_polynom):
            return polynom_multiplication
        else:
            diff = len(polynom_multiplication) - len(irreducible_polynom)
            if diff == 0:
                polynom_multiplication = np.polyadd(
                    polynom_multiplication, irreducible_polynom)
                polynom_multiplication = np.array(
                    [i % 2 for i in polynom_multiplication])  # mod 2 für jede Stelle
                return polynom_multiplication
            else:
                sx_polynom = generate_sx_polynom(diff)
                sx_irr_polynom = np.polymul(sx_polynom, irreducible_polynom)
                # mod 2 für jede Stelle
                sx_irr_polynom = np.array([i % 2 for i in sx_irr_polynom])

                polynom_multiplication = np.polyadd(
                    polynom_multiplication, sx_irr_polynom)
                polynom_multiplication = np.array(
                    [i % 2 for i in polynom_multiplication])  # mod 2 für jede Stelle
                polynom_multiplication = np.trim_zeros(
                    polynom_multiplication, trim='f')


def calculate_multiplication(polynom1, polynom2, irreducible_polynom, e):
    polynom_multiplication = np.polymul(polynom1, polynom2)
    polynom_multiplication = np.array(
        [i % 2 for i in polynom_multiplication])  # mod 2 für jede Stelle
    return calculate_mod(polynom_multiplication, irreducible_polynom, e)


def calculate_multiplication_table(irreducible_polynom):
    e = len(irreducible_polynom)-1
    irreducible_polynom = np.array(
        [int(d) for d in irreducible_polynom.zfill(e)])

    table_string = "x\t"
    for i in range(2**e - 1):
        table_string += str((i+1))+"\t"
    print(table_string)
    print()
    for i in range(1, 2**e):
        print(i, "\t", end='')
        polynom_i = generate_polynom(i, e)
        for j in range(1, 2**e):
            polynom_j = generate_polynom(j, e)
            result = calculate_multiplication(
                polynom_i, polynom_j, irreducible_polynom, e)
            dez = int(str(result)[1:-1].replace(' ', ''), 2)
            print(dez, '\t', end='')
        print()


if __name__ == '__main__':
    # input_polynom = input("Polynom: ")
    input_polynom = "1101"
    calculate_multiplication_table(input_polynom)
