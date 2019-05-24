import unittest
from poplars.hypermut import rate_ratio


class TestRateRatio(unittest.TestCase):
    def testRateRatio(self):
        ctable = [[4, 17], [7, 33]]
        expected = 1.1
        output = rate_ratio(ctable)
        result = round(output, 1)

        self.assertEqual(expected, result)

    def test1UndefinedRateRatio(self):
        ctable = [[4, 17], [0, 33]]
        expected = "undef"
        result = rate_ratio(ctable)
        self.assertEqual(expected, result)

    def test2UndefinedRateRatio(self):
        ctable = [[4, 17], [7, 0]]
        expected = "undef"
        result = rate_ratio(ctable)
        self.assertEqual(expected, result)
