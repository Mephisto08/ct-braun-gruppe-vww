import time


class P:  # Polynomial
    def __init__(self, value):
        self.value = str(value)

    def __add__(p1, p2):
        width = max(len(p1.value), len(p2.value))
        p1, p2 = p1.pad(width), p2.pad(width)
        result = ''
        for i in range(width):
            pv1, pv2 = int(p1[i]), int(p2[i])
            result += str(pv1 + pv2)
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

    def __mod__(p, mod):
        result = ''
        for i in range(len(p.value)):
            pv = int(p.value[i])
            result += str(pv % mod)
        return P(result)

    def reduce(input, irreducible_p, p):
        result = input % p
        diff = len(result.value) - len(irreducible_p.value)
        if diff >= 0:
            while True:
                sx = P('1' + '0' * diff)
                ir_sx = irreducible_p * sx
                result += ir_sx
                result %= p
                result.value = result.value.lstrip('0')  # Remove leading zero, since degree got reduced
                diff = len(result.value) - len(irreducible_p.value)
                if diff < 0:
                    break
        return result

    def pad(p, width) -> str:
        return p.value.zfill(width)


class MulTab:  # Multiplication Table
    def __init__(self, p, irreducible_p):
        self.p = p
        self.e = len(irreducible_p.value) - 1
        self.irreducible_p = irreducible_p
        self.values = self.calc_table()

    def calc_table(self):
        width = self.p ** self.e
        result = [[P(self.bin(0))] * width for w in range(width)]  # Initialize two-dimensional array
        for i in range(1, width):
            for j in range(i, width):
                res = self.mul_mod(P(self.bin(i)), P(self.bin(j)))
                result[i][j] = res
                result[j][i] = res
        return result

    def mul_mod(self, p1, p2):
        mul_p = p1 * p2
        result = mul_p.reduce(self.irreducible_p, self.p)
        return result

    def print(self):
        for row in self.values:
            for field in row:
                # print(self.pad(field.value) + ' ', end='')    # Binary
                print(str(int(field.value, 2)) + ' ', end='')   # Decimal
            print()

    def bin(self, number):
        return self.pad(bin(number)[2:])

    def pad(self, number):
        return number.zfill(self.e)


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
    p = 2
    e = 8

    start_time = time.time()
    mt = MulTab(p, P(ips[e]))
    stop_time = time.time()

    mt.print()
    print("--- %s seconds ---" % (stop_time - start_time))

