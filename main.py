import pandas as pd
import time
import math

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)

ips = [  # Irreducible Polynomial Lookup Table
    '',
    '',
    111,  # 2
    1101,  # 3
    11001,  # 4
    100101,  # 5
    1100001,  # 6
    11000001,  # 7
    100011101,  # 8
    1000010001,  # 9
    10000001001,  # 10
]


class P:  # Polynomial
    def __init__(self, value):
        self.value = str(value)

    def __eq__(p1, p2):
        p1 = p1.value.lstrip('0')  # Remove leading zeros
        p2 = p2.value.lstrip('0')
        return p1 == p2

    def __add__(p1, p2):
        width = max(len(p1.value), len(p2.value))
        p1, p2 = p1.pad(width), p2.pad(width)
        result = ''
        for i in range(width):
            pv1, pv2 = int(p1[i]), int(p2[i])
            result += str(pv1 + pv2)
        return P(result)

    def __sub__(p1, p2):
        width = max(len(p1.value), len(p2.value))
        p1, p2 = p1.pad(width), p2.pad(width)
        result = ''
        for i in range(width):
            pv1, pv2 = int(p1[i]), int(p2[i])
            result += str(pv1 - pv2)
        result = result.lstrip('0')  # Remove leading zeros
        return P(result)

    def __mul__(p1, p2):
        pl1, pl2 = len(p1.value), len(p2.value)
        width = pl1 + pl2 - 1
        result = [0] * width
        for i in range(pl1):
            for j in range(pl2):
                pv1, pv2 = int(p1.value[i]), int(p2.value[j])
                result[i + j] += pv1 * pv2
        result = ''.join(str(x) for x in result)  # List to string
        result = result.lstrip('0')  # Remove leading zeros
        return P(result or '0')

    def __truediv__(p1, p2):
        diff = len(p1.value) - len(p2.value)
        if diff < 0:
            return [P('0'), p1]

        pv1, pv2 = int(p1.value[0]), int(p2.value[0])
        factor = int(pv1 / pv2)
        pr = P(str(factor) + '0' * diff)
        pmul = pr * p2
        psub = (p1 - pmul).abs()

        if psub == P('0'):
            return [pr, P('0')]

        pdiv = psub / p2
        return [pr + pdiv[0], pdiv[1]]

    def __floordiv__(p1, p2):
        return (p1 / p2)[0]

    def __mod__(p1, p2):
        mod = (p1 / p2)[1]
        return mod

    def mul(p, factor: int):
        result = ''
        for i in range(len(p.value)):
            pv = int(p.value[i])
            result += str(int(pv * factor))
        return P(result)

    def div(p, factor: int):
        result = ''
        for i in range(len(p.value)):
            pv = int(p.value[i])
            result += str(int(pv / factor))
        return P(result)

    def mod(p, mod: int):
        result = ''
        for i in range(len(p.value)):
            pv = int(p.value[i])
            result += str(pv % mod)
        return P(result)

    def reduce(input, irreducible_p, p):
        result = input.mod(p)
        diff = len(result.value) - len(irreducible_p.value)
        if diff >= 0:
            while True:
                sx = P('1' + '0' * diff)
                ir_sx = irreducible_p * sx
                result += ir_sx
                result = result.mod(p)
                # Remove leading zero, since degree got reduced
                result.value = result.value.lstrip('0')
                diff = len(result.value) - len(irreducible_p.value)
                if diff < 0:
                    break
        return result

    def pad(p, width) -> str:
        return p.value.zfill(width)

    def abs(self):
        pv = self.value.replace('-', '')
        return P(pv)


class MulTab:  # Multiplication Table
    def __init__(self, irreducible_p, p=2):
        self.p = p
        self.e = len(irreducible_p.value) - 1
        self.irreducible_p = irreducible_p
        if self.e < 2:
            raise ValueError("e cannot be less than 2")
        else:
            self.width = self.p ** self.e
            self.values = [[P(self.to_base(0))] * self.width for w in
                           range(self.width)]  # Initialize two-dimensional array

    def calc_table(self):
        for i in range(1, self.width):
            for j in range(i, self.width):
                res = self.mul_mod(P(self.to_base(i)), P(self.to_base(j)))
                self.values[i][j] = res
                self.values[j][i] = res

    def mul_mod(self, p1, p2):
        mul_p = p1 * p2
        result = mul_p.reduce(self.irreducible_p, self.p)
        return result

    def print(self, raw=False):
        df = pd.DataFrame(self.values)
        def format(field): return int(
            self.pad(field.value), 2)  # Decimal converted
        if raw:
            def format(field): return self.pad(field.value)  # Raw
        df = df.applymap(format)
        print(df)

    def to_base(self, n):
        base = self.p
        digits = ""
        while n:
            digits = str(int(n % base)) + digits
            n //= base  # Floor division
        return self.pad(digits)

    def bin(self, number):
        return self.pad(bin(number)[2:])

    def pad(self, number):
        return number.zfill(self.e)


# Extended Euclidean Algorithm
def eea(p1, p2, irreducible_p, p):
    if p1 == P("0"):
        return p2, P("0"), P("1")
    gcd, u, v = eea(p2 % p1, p1, irreducible_p, p)
    pfdivmul = ((p2 // p1) * u).reduce(irreducible_p, p)
    x = (v - pfdivmul).abs()
    y = u
    return gcd, x, y


def calc_kgv(a, b):
    """ calculate the least common multiple """
    return abs(a * b) // math.gcd(a, b)


def print_matrix(m):
    for p in m:
        print(p.value)


def gen_em(rows: int):
    result = []
    for i in range(rows):
        result.append(bin(2**(rows-i-1))[2:].zfill(rows))
    return result


def gen_transposed_matrix(m):
    t_matrix = []
    for i in range(len(m[0].value)):
        pol_str = ""
        for j in range(len(m)):
            pol_str += m[j].value[i]

        t_matrix.append(P(pol_str))

    return t_matrix


def generate_reducedRowEchelonForm(M, e):
    if not M:
        return
    lead = 0
    rowCount = len(M)
    columnCount = len(M[0].value)
    for r in range(rowCount):
        if lead >= columnCount:
            return M
        i = r

        while M[i].value[lead] == "0":
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return M
        M[i], M[r] = M[r], M[i]
        lv = M[r].value[lead]
        M[r] = P(''.join([str(abs(int(int(mrx) / float(lv))))
                          for mrx in M[r].value])).mod(e)  # !VORSICHT
        for i in range(rowCount):
            if i != r:
                lv = M[i].value[lead]
                M[i] = P(''.join([str(abs(int(iv) - int(lv) * int(rv)))
                                  for rv, iv in zip(M[r].value, M[i].value)])).mod(e)
        lead += 1
    return M


def generate_control_matrix(gm):
    g = gm.copy()
    rowCount = len(g)

    for i in range(rowCount):
        g[i] = P(g[i].value[rowCount:])

    p_transposed = gen_transposed_matrix(g)

    em = gen_em(len(p_transposed))

    h = []
    for i in range(len(p_transposed)):
        h.append(P(p_transposed[i].value + em[i]))

    km = gen_transposed_matrix(h)

    return km


def generate_syndrom_table(km):
    n = len(km)
    syndrom_table = {}

    syndrom_table['0'*len(km[0].value)] = P('0'*n)
    for i in range(n):
        cur_pol = str(bin(2**(i))[2:].zfill(n))
        syndrom_table[km[len(km)-1-i].value] = P(cur_pol)

    for i in range(2**n):
        cur_pol = str(bin(i)[2:].zfill(n))
        temp_pol = P('0'*len(km[0].value))
        for j in range(n):
            if cur_pol[j] == '0':
                continue

            temp_pol += km[j]

        temp_pol = temp_pol.mod(2)

        if temp_pol.value not in syndrom_table:
            syndrom_table[temp_pol.value] = P(cur_pol)

    return syndrom_table


def error_correction_with_syndrom_table(code_polynom, km, syndrom_table):
    n = len(km)

    syndrom_class = P('0'*len(km[0].value))
    for j in range(n):
        if code_polynom.value[j] == '0':
            continue

        syndrom_class += km[j]

    syndrom_class = syndrom_class.mod(2).value
    error_polynom = syndrom_table[syndrom_class]

    return syndrom_class, (code_polynom + error_polynom).mod(2)


def calc_g_mul_ht(gm, km):
    n = len(gm[0].value)

    temp = P('0' * len(km[0].value))
    for e_gm in gm:
        for j in range(n):
            t = e_gm.value[j]
            if e_gm.value[j] == '0':
                continue

            temp += km[j]

    result = temp.mod(2).value

    return P(result)


# Aufgabe 4
def is_power_of_two(n):
    return (n != 0) and (n & (n-1) == 0)


def generate_hamming_control_matrix(m):
    control_matrix = []
    for i in range(1, 2**m):
        if not is_power_of_two(i):
            item = P(str(bin(i))[2:].zfill(m))
            control_matrix.append(item)

    for i in gen_em(m):
        control_matrix.append(P(i))

    return control_matrix


def hamming_control_matrix_to_generator_matrix(control_matrix):
    m = len(control_matrix[0].value)

    generator_matrix = []

    for i in gen_em(len(control_matrix)-m):
        generator_matrix.append(P(i))

    p_transposed = gen_transposed_matrix(
        control_matrix[:-m])

    generator_matrix += p_transposed

    return generator_matrix


def decode_hamming(codeword, control_matrix):
    cm_transposed = gen_transposed_matrix(control_matrix)

    # y * H^T berechnen
    pol_str = ""
    for i in cm_transposed:
        sum = 0
        for j in range(len(codeword.value)):
            sum += int(codeword.value[j]) * int(i.value[j])
        pol_str += str(sum)

    y_ht = P(pol_str).mod(2)

    # faktor a berechnen
    a = 0
    for i in y_ht.value:
        if int(i) > 0:
            a = int(i)
            break

    # fehler berechnen
    error = None
    for i, item in enumerate(gen_transposed_matrix(cm_transposed)):
        if item.mul(a) == y_ht:
            error_string = '0'*len(cm_transposed[0].value)
            error_string = error_string[:i] + str(a) + error_string[i+1:]
            error = P(error_string)

    if error:
        corrected_codeword = (codeword + error).mod(2)
        return corrected_codeword

    return codeword


if __name__ == '__main__':
    # Choose an e between 2 and 8
    e = 4

    print("\n┎────────────────┒")
    print("┃   Exercise 1   ┃")
    print("┖────────────────┚")
    print("e = " + str(e))
    start_time = time.time()
    mt = MulTab(P(ips[e]))
    mt.calc_table()
    stop_time = time.time()
    mt.print()
    print("\n➜ Took %s seconds\n" % (stop_time - start_time))

    print("\n┎────────────────┒")
    print("┃   Exercise 2   ┃")
    print("┖────────────────┚")
    print("e = " + str(e))
    start_time = time.time()
    df = pd.DataFrame()
    df.index = ['Field element', 'GDC', 'u', 'v', 'mul result']
    for i in range(1, 2 ** e):
        p1 = P(bin(i)[2:])
        gcd, u, v = eea(p1, mt.irreducible_p, mt.irreducible_p, mt.p)
        mul_r = mt.mul_mod(p1, u)
        result = [p1.value, gcd.value, u.value, v.value, mul_r.value]
        df = df.assign(**{str(i): result})
    stop_time = time.time()
    print(df)
    print("\n➜ Took %s seconds\n" % (stop_time - start_time))

    print("\n┎────────────────┒")
    print("┃   Exercise 3   ┃")
    print("┖────────────────┚")
    print("e = " + str(e))
    start_time = time.time()

    gm = [
        P("11010"),
        P("11010")
    ]

    print("\nGenerator-Matrix:")
    print_matrix(gm)
    print()

    kgm = generate_reducedRowEchelonForm(gm, 2)

    print("\nKanonische-Generator-Matrix:")
    print_matrix(kgm)
    print()

    km = generate_control_matrix(kgm)

    print("\nKontroll-Matrix:")
    print_matrix(km)
    print()

    syndrom_table = generate_syndrom_table(km)

    print("\nSyndrom Tabelle:\nSyndr.\tError")
    for key, value in syndrom_table.items():
        print(key, '\t', value.value)
    print()

    corrected_codeword = error_correction_with_syndrom_table(
        P("01101"), km, syndrom_table)
    print("Syndrom class:", corrected_codeword[0])
    print("Corrected codeword with syndrom table:",
          corrected_codeword[1].value)

    g_mul_ht_result = calc_g_mul_ht(gm, km)
    print("G * Ht: ", g_mul_ht_result.value)

    stop_time = time.time()
    print("\n➜ Took %s seconds\n" % (stop_time - start_time))

    print("\n┎────────────────┒")
    print("┃   Exercise 4   ┃")
    print("┖────────────────┚")

    # Choose an m greater or equal 3
    m = 3

    km = generate_hamming_control_matrix(m)

    print("Kontroll-Matrix:")
    for i in range(len(km[0].value)):
        for j in km:
            print(j.value[i], end=' ')
        print()

    gm = hamming_control_matrix_to_generator_matrix(km)

    print("\nGenerator-Matrix:")
    for i in range(len(gm[0].value)):
        for j in gm:
            print(j.value[i], end=' ')
        print()

    # Choose codeword containing 0 and 1 with length (2^m - 1)
    codeword = P("0101111")
    corrected_codeword = decode_hamming(codeword, km)
    print("Received Codeword:", codeword.value)
    print("Korrigiertes Codeword:", corrected_codeword.value)
