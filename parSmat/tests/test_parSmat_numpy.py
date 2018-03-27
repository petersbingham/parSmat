import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import parSmat as psm
import channelutil as chanutil

import unittest

class parentTest(unittest.TestCase):
    def calculateCoefficients(self, thres, sMatData):
        psm.usePythonTypes()
        asymCal = chanutil.asymCal(chanutil.HARTs, thresholds=thres)
        return psm.calculateCoefficients(sMatData, asymCal), asymCal

class test_parSmat(parentTest):
    def runTest(self):
        import numpyTestData as dat
        coeffs,_ = self.calculateCoefficients([0.,2.],dat.sMatData_inel)
        for i in range(4):
            testdps = 1e-9
            self.assertTrue(psm.nw.areMatricesClose(dat.A_inel[i],coeffs[0][i],
                                                    rtol=testdps, atol=testdps))
            self.assertTrue(psm.nw.areMatricesClose(dat.B_inel[i],coeffs[1][i],
                                                    rtol=testdps, atol=testdps))

class test_fin(parentTest):
    def runTest(self):
        import numpyTestData as dat
        coeffs, asymCal = self.calculateCoefficients([0.,0.],dat.sMatData_el)
        fun = psm.getElasticFinFun(coeffs, asymCal)
        parFinMat = fun(3.0)
        testdps = 1e-9
        self.assertTrue(psm.nw.areMatricesClose(parFinMat,dat.finData_el_3,
                                                rtol=testdps, atol=testdps))

class test_Smat(parentTest):
    def runTest(self):
        import numpyTestData as dat
        coeffs, asymCal = self.calculateCoefficients([0.,0.],dat.sMatData_el)
        fun = psm.getElasticSmatFun(coeffs, asymCal)
        parSmat = fun(3.0)
        testdps = 1e-9
        self.assertTrue(psm.nw.areMatricesClose(parSmat,dat.sMatData_el_3,
                                                rtol=testdps, atol=testdps))
        

if __name__ == "__main__":
    #Just for debug
    b = test_Smat()
    b.runTest()
