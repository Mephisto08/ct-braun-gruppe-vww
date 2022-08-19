import pandas as pd
import time

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
        if self.e < 2:
            raise ValueError("e cannot be less than 2")
        else:
            self.width = self.p ** self.e
            self.values = [[P(self.to_base(0))] * self.width for w in range(self.width)]  # Initialize two-dimensional array

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
        format = lambda field: int(self.pad(field.value), 2)  # Decimal converted
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
def eea(p1, p2, irreducible_p, p):
    if p1 == P("0"):
        return p2, P("0"), P("1")
    gcd, u, v = eea(p2 % p1, p1, irreducible_p, p)
    pfdivmul = ((p2 // p1) * u).reduce(irreducible_p, p)
    x = (v - pfdivmul).abs()
    y = u
    return gcd, x, y


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
    for i in range(1, 2**e):
        p1 = P(bin(i)[2:])
        gcd, u, v = eea(p1, mt.irreducible_p, mt.irreducible_p, mt.p)
        mul_r = mt.mul_mod(p1, u)
        result = [p1.value, gcd.value, u.value, v.value, mul_r.value]
        df = df.assign( **{str(i): result})
    stop_time = time.time()
    print(df)
    print("\n➜ Took %s seconds\n" % (stop_time - start_time))
