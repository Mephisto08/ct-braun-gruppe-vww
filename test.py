import unittest
from main import *

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


class MyTestCase(unittest.TestCase):
    def test_prim_2_poly_3(self):
        p = 2
        e = 3
        testMulTab_3 = MulTab(p, P(ips[e]))
        # test position 7x7
        self.assertEqual(testMulTab_3.values[7][7].value, "10")

        # test position 6x4
        self.assertEqual(testMulTab_3.values[6][4].value, "111")

        # check if mirroring of the matrix works correctly
        self.assertEqual(testMulTab_3.values[4][6].value, "10")


if __name__ == '__main__':
    unittest.main()
