import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import parSmat as psm
import channelutil as chanutil

import unittest

class test_parSmat(unittest.TestCase):
    def runTest(self):
        import numpyCoeffsData as dat
        psm.usePythonTypes()
        chanCalc = chanutil.asymCal([0.,2.],units=chanutil.hartrees)
        coeffs = psm.calculateCoefficients(dat.sMatData, chanCalc)
        for i in range(4):
            testdps = 1e-9
            self.assertTrue(psm.nw.areMatricesClose(dat.A[i],coeffs[0][i],rtol=testdps, atol=testdps))
            self.assertTrue(psm.nw.areMatricesClose(dat.B[i],coeffs[1][i],rtol=testdps, atol=testdps))

if __name__ == "__main__":
    #Just for debug
    b = test_parSmat()
    b.runTest()
