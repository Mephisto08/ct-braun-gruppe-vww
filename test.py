import unittest
from main import *


def test_eea_e(test, e):
    mt = MulTab(P(ips[e]))
    for i in range(1, 2 ** e):
        p1 = P(bin(i)[2:])
        gcd, u, v = eea(p1, mt.irreducible_p, mt.irreducible_p, mt.p)
        mul_r = mt.mul_mod(p1, u)
        test.assertEqual("1", mul_r.value)

        # Check for missing modulo
        values = [p1, gcd, u, v, mul_r]
        for v in values:
            for c in v.value:
                test.assertLess(int(c), int(mt.p))


def test_multab_e(test, e):
    mt = MulTab(P(ips[e]))
    mt.calc_table()

    for i in range(1, mt.width):
        row_one_counter = 0
        column_one_counter = 0
        for j in range(1, mt.width):
            test.assertEqual(mt.values[i][j].value, mt.values[j][i].value)  # Check symmetry of matrix

            # Check for missing modulo
            field = mt.values[i][j].value
            for c in field:
                test.assertLess(int(c), int(mt.p))

            v_row = mt.values[i][j].value.lstrip('0')  # Remove leading zero
            v_column = mt.values[j][i].value.lstrip('0')
            if v_row == '1':
                row_one_counter += 1
            if v_column == '1':
                column_one_counter += 1
        test.assertEqual(1, row_one_counter)  # Check count of found one's in row
        test.assertEqual(1, column_one_counter)  # Check count of found one's in column


class TestCase(unittest.TestCase):

    def test_multab(self):
        for e in range(2, 9):
            print("Running test_multab with e =", e)
            test_multab_e(self, e)

    def test_eea(self):
        for e in range(2, 9):
            print("Running test_eea with e =", e)
            test_eea_e(self, e)


if __name__ == '__main__':
    unittest.main()
