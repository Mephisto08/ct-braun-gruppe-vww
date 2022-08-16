import unittest
from main import *


class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(hi('Tester'), "Hi, Tester!!!!!!!!!")


if __name__ == '__main__':
    unittest.main()
