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
        for c in p.value:
            print(c + ' ', end='')
        print()


def print_dec_matrix(m):
    for p in m:
        for c in p:
            print(c + ' ', end='')
        print()


def dec_to_bin_pol(e, dec):
    binary = ""
    for d in dec:
        b = bin(int(d))[2:].zfill(e)
        binary += b
    return P(binary)


def dec_array_to_bin_pol_array(e, dec_array):
    bin_pol_array = []
    for d in dec_array:
        p = dec_to_bin_pol(e, d)
        bin_pol_array.append(p)
    return bin_pol_array


def bin_pol_to_dec(e, bin_pol):
    dec = ""
    for i in range(0, len(bin_pol.value), e):
        dec += str(int(bin_pol.value[i:i+e], 2))
    return dec


def bin_pol_array_to_dec_array(e, bin_pol_array):
    dec_array = []
    for p in bin_pol_array:
        d = bin_pol_to_dec(e, p)
        dec_array.append(d)
    return dec_array


def gen_em(rows: int):
    result = []
    for i in range(rows):
        result.append(bin(2**(rows-i-1))[2:].zfill(rows))
    return result


def gen_transposed_matrix(m):
    t_matrix = []
    for i in range(len(m[0])):
        pol_str = ""
        for j in range(len(m)):
            pol_str += m[j][i]

        t_matrix.append(pol_str)

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


def generate_canonical_generator_matrix(M, e):
    rref = generate_reducedRowEchelonForm(M, e)
    result = []
    for i in rref:
        if '1' in i.value:
            result.append(i)
    return result


def generate_control_matrix(gm):
    g = gm.copy()
    rowCount = len(g)

    for i in range(rowCount):
        g[i] = g[i][rowCount:]

    p_transposed = gen_transposed_matrix(g)

    em = gen_em(len(p_transposed))

    h = []
    for i in range(len(p_transposed)):
        h.append(p_transposed[i] + em[i])

    km = gen_transposed_matrix(h)

    return km


def generate_syndrom_table(km, q):
    n = len(km)
    syndrom_table = {}

    syndrom_table['0'*len(km[0])] = '0'*n
    for i in range(n):
        for j in range(q):
            cur_pol = str(bin(2**(i))[2:].zfill(n))
            syndrom_table[km[len(km)-1-i]] = P(cur_pol)

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


# Aufgabe 5
def generate_reed_muller_code(r, m):
    if r == 0:
        return [P('1' * 2**m)]
    elif r > m:
        return generate_reed_muller_code(m, m)

    rm_1 = generate_reed_muller_code(r, m-1)
    rm_2 = generate_reed_muller_code(r-1, m-1)

    rm_generator_matrix = []

    for i in range(len(rm_1)):
        rm_generator_matrix.append(P(rm_1[i].value + rm_1[i].value))

    for i in range(len(rm_2)):
        rm_generator_matrix.append(P('0' * len(rm_1[0].value) + rm_2[i].value))

    return rm_generator_matrix


# Aufgabe 6
def determine_primitive_element(q):
    gf_target = [x for x in range(1, q)]

    for alpha in range(1, q):
        gf_without_zero = []
        for i in range(q-1): # 0 <= i <= q-2
            gf_without_zero.append((alpha ** i) % q)

        if set(gf_without_zero) == set(gf_target):
            return alpha

    return None


def add_with_mod(p1, p2, q):
    width = max(len(p1.value), len(p2.value))
    p1, p2 = p1.pad(width), p2.pad(width)
    result = ''
    for i in range(width):
        pv1, pv2 = int(p1[i]), int(p2[i])
        result += str((pv1 + pv2) % q)
    return P(result)


def mul_with_mod(p1, p2, q):
    pl1, pl2 = len(p1.value), len(p2.value)
    width = pl1 + pl2 - 1
    result = [0] * width
    for i in range(pl1):
        for j in range(pl2):
            pv1, pv2 = int(p1.value[i]), int(p2.value[j])
            result[i + j] = (result[i + j] + (pv1 * pv2) % q) % q
    result = ''.join(str(x) for x in result)  # List to string
    result = result.lstrip('0')  # Remove leading zeros
    return P(result or '0')


def generate_reed_solomon_generator_polynom(alpha, q, d):
    g_list = []
    for i in range(1, d): # 1 <= i <= d-1
        value = -(alpha ** i) % q
        p = P('1' + str(value))
        g_list.append(p)

    result_polynom = P('1')
    for i in range(len(g_list)):
        result_polynom = mul_with_mod(result_polynom, g_list[i], q)

    return result_polynom


def generate_reed_solomon_control_polynom(alpha, q, d):
    g_list = [P('1' + str(-1 % q))]
    for i in range(d, q-1): # d <= i <= q-2
        value = -(alpha ** i) % q
        p = P('1' + str(value))
        g_list.append(p)

    result_polynom = P('1')
    for i in range(len(g_list)):
        result_polynom = mul_with_mod(result_polynom, g_list[i], q)

    return result_polynom


def generate_reed_solomon_control_matrix(control_polynom, d):
    H = []
    max = d-1
    for i in range(d-1):
        pol_str = ('0' * i) + control_polynom.value + ('0' * (max-i-1))
        H.append(P(pol_str))
    return H


def generate_reed_solomon_generator_matrix(generator_polynom, q, d):
    inverted_gp_string = generator_polynom.value[::-1]
    n = q-1

    G = []
    max = n - len(inverted_gp_string) + 1
    for i in range(max):
        pol_str = ('0' * i) + inverted_gp_string + ('0' * (max-i-1))
        G.append(P(pol_str))
    return G


def generate_reed_solomon_vandermonde_matrix(polynom, q):
    n = len(polynom.value)
    m = n

    V = []
    for i in range(0, m):  # 0 <= i <= m-1
        p_string = ''
        for j in range(n):  # 0 <= j <= n-1
            value = (int(polynom.value[i]) ** j) % q
            p_string += str(value)
        V.append(P(p_string))

    return V


def generate_reed_solomon_code(e, d):
    alpha = determine_primitive_element(2**e)

    #generator_polynom = generate_reed_solomon_generator_polynom(alpha, q, d)
    #generator_matrix = generate_reed_solomon_generator_matrix(generator_polynom, q, d)
    #control_polynom = generate_reed_solomon_control_polynom(alpha, q, d)
    #control_matrix = generate_reed_solomon_control_matrix(control_polynom, d)
    #vandermonde_matrix = generate_reed_solomon_vandermonde_matrix(generator_polynom, q)

    print("e:", e)
    print("d:", d)

    #print("\nAlpha:", alpha)
    #print("Generator-Polynom:", generator_polynom.value)
    #print("Kontroll-Polynom:", control_polynom.value)

    #print("\nGenerator-Matrix:")
    #print_matrix(generator_matrix)

    #print("\nKontroll-Matrix:")
    #print_matrix(control_matrix)

    #print("\nVandermonde-Matrix:")
    #print_matrix(vandermonde_matrix)


def exercise1():
    print("\n┎────────────────┒")
    print("┃   Exercise 1   ┃")
    print("┖────────────────┚")
    print("» Multiplikationstabelle «")

    # Choose an e between 2 and 8
    e = 4

    start_time = time.time()
    mt = MulTab(P(ips[e]))
    mt.calc_table()
    stop_time = time.time()

    print("e: ", e, "\n")
    mt.print()
    print("\n➜ In %s Sekunden" % (stop_time - start_time))


def exercise2():
    print("\n┎────────────────┒")
    print("┃   Exercise 2   ┃")
    print("┖────────────────┚")
    print("» Erweiterter Euklidischer Algorithmus «")

    # Choose an e between 2 and 8
    e = 4

    start_time = time.time()
    mt = MulTab(P(ips[e]))
    df = pd.DataFrame()
    df.index = ['Field element', 'GDC', 'u', 'v', 'mul result']

    for i in range(1, 2 ** e):
        p1 = P(bin(i)[2:])
        gcd, u, v = eea(p1, mt.irreducible_p, mt.irreducible_p, mt.p)
        mul_r = mt.mul_mod(p1, u)
        result = [p1.value, gcd.value, u.value, v.value, mul_r.value]
        df = df.assign(**{str(i): result})

    stop_time = time.time()

    print("e: ", e, "\n")
    print(df)
    print("\n➜ In %s Sekunden" % (stop_time - start_time))


def exercise3():
    print("\n┎────────────────┒")
    print("┃   Exercise 3   ┃")
    print("┖────────────────┚")
    print("» Linearer-Code «\n")

    # Choose an e between 2 and 8
    e = 2

    # Choose generator matrix
    gm = [
        "10111",
        "01123"
    ]

    # Choose received codeword
    codeword = "11110"

    n = len(gm[0])
    print("GF(2^" + str(e) + ")^" + str(n) + " = GF(" + str(2 ** e) + ")^" + str(n))
    print("\nGenerator-Matrix:")
    print_dec_matrix(gm)
    print("\nGenerator-Matrix (Binär):")
    print_matrix(dec_array_to_bin_pol_array(e, gm))

    gm = dec_array_to_bin_pol_array(e, gm)
    codeword = dec_to_bin_pol(e, codeword)

    # kanonische
    kgm = generate_canonical_generator_matrix(gm, 2)
    dec_kgm = bin_pol_array_to_dec_array(e, kgm)

    # kontroll
    km = generate_control_matrix(dec_kgm)
    dec_km = dec_array_to_bin_pol_array(e, km)

    syndrom_table = generate_syndrom_table(km, 2**e)
    syndrom_class, corrected_codeword = error_correction_with_syndrom_table(codeword, km, syndrom_table)
    g_mul_ht_result = calc_g_mul_ht(gm, km)

    print("\nKanonische-Generator-Matrix:")
    print_dec_matrix(bin_pol_array_to_dec_array(e, kgm))
    print("\nKanonische-Generator-Matrix (Binär):")
    print_matrix(kgm)

    print("\nKontroll-Matrix:")
    print_dec_matrix(bin_pol_array_to_dec_array(e, km))
    print("\nKontroll-Matrix (Binär):")
    print_matrix(km)

    print("\nSyndrom Tabelle:\nSyndr.\tError")
    for key, value in syndrom_table.items():
        print(key, '\t', value.value, sep='')

    print("\nSyndrom Klasse:", syndrom_class)
    print("Empfangenes Codeword:", bin_pol_to_dec(e, codeword), "(Binär: " + str(codeword.value) + ")")
    print("Korrigiertes Codeword:", bin_pol_to_dec(e, corrected_codeword), "(Binär: " + str(corrected_codeword.value) + ")")
    print("G * Ht:", bin_pol_to_dec(e, g_mul_ht_result), "(Binär: " + str(g_mul_ht_result.value) + ")")


def exercise4():
    print("\n┎────────────────┒")
    print("┃   Exercise 4   ┃")
    print("┖────────────────┚")
    print("» Hamming-Code «")

    # Choose an m greater or equal 3
    m = 3

    # Choose received codeword containing 0 and 1 with length (2^m - 1)
    codeword = P("0101111")

    km = generate_hamming_control_matrix(m)
    gm = hamming_control_matrix_to_generator_matrix(km)
    corrected_codeword = decode_hamming(codeword, km)

    print("m:", m)

    print("\nKontroll-Matrix:")
    for i in range(len(km[0].value)):
        for j in km:
            print(j.value[i], end=' ')
        print()

    print("\nGenerator-Matrix:")
    for i in range(len(gm[0].value)):
        for j in gm:
            print(j.value[i], end=' ')
        print()

    print("\nEmpfangenes Codeword:", codeword.value)
    print("Korrigiertes Codeword:", corrected_codeword.value)


def exercise5():
    print("\n┎────────────────┒")
    print("┃   Exercise 5   ┃")
    print("┖────────────────┚")
    print("» Reed-Muller-Code «")

    # Choose r, m for Reed-Muller-Code construction
    r = 1
    m = 5

    reed_muller_code = generate_reed_muller_code(r, m)

    print("r:", r)
    print("m:", m)

    print("\nGenerator-Matrix:")
    print_matrix(reed_muller_code)


def exercise6():
    print("\n┎────────────────┒")
    print("┃   Exercise 6   ┃")
    print("┖────────────────┚")
    print("» Reed-Solomon-Code «")

    # Choose e, d for Reed-Solomon-Code construction with q = 2^e
    e = 2
    d = 5

    generate_reed_solomon_code(e, d)


if __name__ == '__main__':
    #exercise1()
    #exercise2()
    exercise3()
    #exercise4()
    #exercise5()
    #exercise6()