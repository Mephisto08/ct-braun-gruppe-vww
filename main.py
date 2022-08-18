import pandas as pd
import time

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)


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
                result.value = result.value.lstrip('0')  # Remove leading zero, since degree got reduced
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
        self.values = self.calc_table()

    def calc_table(self):
        width = self.p ** self.e
        result = [[P(self.to_base(0))] * width for w in range(width)]  # Initialize two-dimensional array
        for i in range(1, width):
            for j in range(i, width):
                res = self.mul_mod(P(self.to_base(i)), P(self.to_base(j)))
                result[i][j] = res
                result[j][i] = res
        return result

    def mul_mod(self, p1, p2):
        mul_p = p1 * p2
        result = mul_p.reduce(self.irreducible_p, self.p)
        return result

    def print(self, raw=False):
        df = pd.DataFrame(self.values)
        format = lambda field: int(field.value, 2)  # Decimal converted
        if raw:
            format = lambda field: self.pad(field.value)  # Raw
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
def eea(p1, p2):
    if p1 == P("0"):
        return p2, P("0"), P("1")

    gcd, u, v = eea(p2 % p1, p1)
    x = (v - (p2 // p1) * u).abs()
    y = u
    return gcd, x, y


if __name__ == '__main__':
    ips = [
        None,
        None,
        110,  # 2
        1101,  # 3
        11001,  # 4
        100101,  # 5
        1100001,  # 6
        11000001,  # 7
        100011101,  # 8
        1000010001,  # 9
        10000001001,  # 10
    ]
    e = 3

    print("###########")
    print("Praktikum 1")
    print("###########\n")
    start_time = time.time()
    mt = MulTab(P(ips[e]))
    stop_time = time.time()
    mt.print()
    print("--- %s seconds ---\n" % (stop_time - start_time))

    print("###########")
    print("Praktikum 2")
    print("###########\n")
    p2 = P(ips[e])
    for i in range(1, 2**e):
        p1 = P(bin(i)[2:])
        gcd, u, v = eea(p1, p2)
        mul_r = mt.mul_mod(p1, u)
        print("Field element:", p1.value, "\n\tGDC:", gcd.value, "\n\tu:\t", u.value, "\n\tv:\t", v.value, "\n\tmul_r:\t", mul_r.value, " (Should be 1)\n")

