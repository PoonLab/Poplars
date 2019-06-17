import unittest

from poplars.riplike import pdistance

class testRiplike(unittest.TestCase):
    def testPdistSimple(self):
        result = pdistance('ACGT', 'ACGC')
        expected = (1, 4)
        self.assertEqual(expected, result)
        
        result = pdistance('ACGT', 'TGCA')
        expected = (4, 4)
        self.assertEqual(expected, result)
    
    def testPdistGapped(self):
        result = pdistance('ACGT', '---T')
        expected = (0, 1)
        self.assertEqual(expected, result)
        
        result = pdistance('ACGT', 'G---')
        expected = (1, 1)
        self.assertEqual(expected, result)
    
    def testBootstrap(self):
        pass
        
