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
            # Check symmetry of matrix
            test.assertEqual(mt.values[i][j].value, mt.values[j][i].value)

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
        # Check count of found one's in row
        test.assertEqual(1, row_one_counter)
        # Check count of found one's in column
        test.assertEqual(1, column_one_counter)


class TestCase(unittest.TestCase):

    # Aufgabe 1
    def test_multab(self):
        print("\n┎────────────────┒")
        print("┃   Exercise 1   ┃")
        print("┖────────────────┚")

        for e in range(2, 9):
            print("Running test_multab with e =", e)
            test_multab_e(self, e)

    # Aufgabe 2
    def test_eea(self):
        print("\n┎────────────────┒")
        print("┃   Exercise 2   ┃")
        print("┖────────────────┚")

        for e in range(2, 9):
            print("Running test_eea with e =", e)
            test_eea_e(self, e)

    # Aufgabe 3
    gm1 = [
        P("11010"),
        P("11010")
    ]

    def test_generate_canonical_generator_matrix(self):
        print("\n┎────────────────┒")
        print("┃   Exercise 3   ┃")
        print("┖────────────────┚")
        print("Running test_generate_canonical_generator_matrix")

        kgm = generate_canonical_generator_matrix(self.gm1, 2)
        self.assertEqual(kgm[0].value, "11010")

    def test_generate_control_matrix(self):
        print("Running test_generate_control_matrix")
        kgm = generate_canonical_generator_matrix(self.gm1, 2)
        km = generate_control_matrix(kgm)
        self.assertEqual(km[0].value, "1010")
        self.assertEqual(km[1].value, "1000")
        self.assertEqual(km[2].value, "0100")
        self.assertEqual(km[3].value, "0010")
        self.assertEqual(km[4].value, "0001")

    def test_generate_syndrom_table(self):
        print("Running test_generate_syndrom_table")
        kgm = generate_canonical_generator_matrix(self.gm1, 2)
        km = generate_control_matrix(kgm)
        syndrom_table = generate_syndrom_table(km)
        self.assertEqual(syndrom_table["0000"].value, "00000")
        self.assertEqual(syndrom_table["0001"].value, "00001")
        self.assertEqual(syndrom_table["0010"].value, "00010")
        self.assertEqual(syndrom_table["0100"].value, "00100")
        self.assertEqual(syndrom_table["1000"].value, "01000")
        self.assertEqual(syndrom_table["1010"].value, "10000")
        self.assertEqual(syndrom_table["0011"].value, "00011")
        self.assertEqual(syndrom_table["0101"].value, "00101")
        self.assertEqual(syndrom_table["0110"].value, "00110")
        self.assertEqual(syndrom_table["0111"].value, "00111")
        self.assertEqual(syndrom_table["1001"].value, "01001")
        self.assertEqual(syndrom_table["1011"].value, "01011")
        self.assertEqual(syndrom_table["1100"].value, "01100")
        self.assertEqual(syndrom_table["1101"].value, "01101")
        self.assertEqual(syndrom_table["1110"].value, "01110")
        self.assertEqual(syndrom_table["1111"].value, "01111")

    def test_error_correction_with_syndrom_table(self):
        print("Running test_error_correction_with_syndrom_table")
        kgm = generate_canonical_generator_matrix(self.gm1, 2)
        km = generate_control_matrix(kgm)
        syndrom_table = generate_syndrom_table(km)
        corrected_codeword = error_correction_with_syndrom_table(
            P("11110"), km, syndrom_table)
        self.assertEqual(corrected_codeword[1].value, "11010")

    def test_g_mul_ht(self):
        print("Running test_g_mul_ht")
        kgm = generate_canonical_generator_matrix(self.gm1, 2)
        km = generate_control_matrix(kgm)
        g_mul_ht_result = calc_g_mul_ht(self.gm1, km)
        self.assertEqual(g_mul_ht_result.value, "0000")

    # Aufgabe 4
    def test_hamming_codes_m3(self):
        print("\n┎────────────────┒")
        print("┃   Exercise 4   ┃")
        print("┖────────────────┚")
        print("Running test_hamming_codes_m3 with m = 3")
        m = 3
        km_m3 = generate_hamming_control_matrix(m)

        codeword_m3 = P("0101111")
        corrected_codeword_m3 = decode_hamming(codeword_m3, km_m3)
        self.assertEqual(corrected_codeword_m3.value, "0001111")

    # Aufgabe 5
    def test_reed_muller_code(self):
        print("\n┎────────────────┒")
        print("┃   Exercise 5   ┃")
        print("┖────────────────┚")

        print("Running test_reed_muller_code with r = 1, m = 3")
        reed_muller_code = generate_reed_muller_code(1, 3)
        self.assertEqual(reed_muller_code[0].value, "11111111")
        self.assertEqual(reed_muller_code[1].value, "01010101")
        self.assertEqual(reed_muller_code[2].value, "00110011")
        self.assertEqual(reed_muller_code[3].value, "00001111")

        print("Running test_reed_muller_code with r = 1, m = 5")
        reed_muller_code = generate_reed_muller_code(1, 5)
        self.assertEqual(
            reed_muller_code[0].value, "11111111111111111111111111111111")
        self.assertEqual(
            reed_muller_code[1].value, "01010101010101010101010101010101")
        self.assertEqual(
            reed_muller_code[2].value, "00110011001100110011001100110011")
        self.assertEqual(
            reed_muller_code[3].value, "00001111000011110000111100001111")
        self.assertEqual(
            reed_muller_code[4].value, "00000000111111110000000011111111")
        self.assertEqual(
            reed_muller_code[5].value, "00000000000000001111111111111111")


if __name__ == '__main__':
    unittest.main()
