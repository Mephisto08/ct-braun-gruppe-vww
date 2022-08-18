import unittest
from main import *


def test_eea_mulr(test, e):
    mt = MulTab(P(ips[e]))
    p2 = P(ips[e])
    for i in range(1, 2 ** e):
        p1 = P(bin(i)[2:])
        gcd, u, v = eea(p1, p2)
        mul_r = mt.mul_mod(p1, u)
        test.assertEqual("1", mul_r.value)


class TestCase(unittest.TestCase):

    def test_prim_2_poly_3(self):
        testMulTab_3 = MulTab(P(ips[3]))
        testMulTab_3.calc_table()
        self.assertEqual("10", testMulTab_3.values[7][7].value)   # test position 7x7
        self.assertEqual("111", testMulTab_3.values[6][3].value)  # test position 6x3
        self.assertEqual("111", testMulTab_3.values[3][6].value)  # check if mirroring of the matrix works correctly

    def test_eea(self):
        for e in range(2, 8):
            test_eea_mulr(self, e)


if __name__ == '__main__':
    unittest.main()
