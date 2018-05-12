import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import parsmat as psm
import channelutil as chanutil

import unittest

class parentTest(unittest.TestCase):
    def calculateCoefficients(self, dat, thres, sMatData):
        psm.useMpmathTypes(dat.TESTDPS)
        chanutil.useMpmathTypes(dat.TESTDPS)
        asymcalc = chanutil.AsymCalc(chanutil.HARTs, thresholds=thres)
        return psm.calculateCoefficients(sMatData, asymcalc), asymcalc

class test_parsmat(parentTest):
    def runTest(self):
        import mpmathTestData as dat
        coeffs,_ = self.calculateCoefficients(dat, [0.,2.], dat.sMatData_inel)
        for i in range(4):
            testdps = 1e-37
            self.assertTrue(psm.nw.areMatricesClose(dat.A_inel[i],coeffs[0][i], 
                                                    rtol=testdps, atol=testdps))
            self.assertTrue(psm.nw.areMatricesClose(dat.B_inel[i],coeffs[1][i], 
                                                    rtol=testdps, atol=testdps))

class test_fin(parentTest):
    def runTest(self):
        import mpmathTestData as dat
        coeffs,asymcalc = self.calculateCoefficients(dat,[0.,0.],dat.sMatData_el)
        fun = psm.getElasticFinFun(coeffs, asymcalc)
        parFinMat = fun(3.0)
        testdps = 1e-37
        self.assertTrue(psm.nw.areMatricesClose(parFinMat,dat.finData_el_3,
                                                rtol=testdps, atol=testdps))

class test_Smat(parentTest):
    def runTest(self):
        import mpmathTestData as dat
        coeffs,asymcalc = self.calculateCoefficients(dat,[0.,0.],dat.sMatData_el)
        fun = psm.getElasticSmatFun(coeffs, asymcalc)
        parsmat = fun(3.0)
        testdps = 1e-37
        self.assertTrue(psm.nw.areMatricesClose(parsmat,dat.sMatData_el_3,
                                                rtol=testdps, atol=testdps))

if __name__ == "__main__":
    #Just for debug
    b = test_parsmat()
    b.runTest()
