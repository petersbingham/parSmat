import os
import sys
basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../..')

import parsmat as psm
import channelutil as chanutil

import unittest

class parentTest(unittest.TestCase):
    def calculate_coefficients(self, thres, smatdata):
        psm.use_python_types()
        asymcalc = chanutil.AsymCalc(chanutil.HARTs, thresholds=thres)
        return psm.calculate_coefficients(smatdata, asymcalc), asymcalc

class test_parsmat(parentTest):
    def runTest(self):
        import numpyTestData as dat
        coeffs,_ = self.calculate_coefficients([0.,2.],dat.smatdata_inel)
        for i in range(4):
            testdps = 1e-9
            self.assertTrue(psm.nw.are_matrices_close(dat.A_inel[i],coeffs[0][i],
                                                    rtol=testdps, atol=testdps))
            self.assertTrue(psm.nw.are_matrices_close(dat.B_inel[i],coeffs[1][i],
                                                    rtol=testdps, atol=testdps))

class test_fin(parentTest):
    def runTest(self):
        import numpyTestData as dat
        coeffs, asymcalc = self.calculate_coefficients([0.,0.],dat.smatdata_el)
        fun = psm.get_elastic_Fin_fun(coeffs, asymcalc)
        parFinMat = fun(3.0)
        testdps = 1e-9
        self.assertTrue(psm.nw.are_matrices_close(parFinMat,dat.finData_el_3,
                                                rtol=testdps, atol=testdps))

class test_Smat(parentTest):
    def runTest(self):
        import numpyTestData as dat
        coeffs, asymcalc = self.calculate_coefficients([0.,0.],dat.smatdata_el)
        fun = psm.get_elastic_Smat_fun(coeffs, asymcalc)
        parsmat = fun(3.0)
        testdps = 1e-9
        self.assertTrue(psm.nw.are_matrices_close(parsmat,dat.smatdata_el_3,
                                                rtol=testdps, atol=testdps))
        

if __name__ == "__main__":
    #Just for debug
    b = test_Smat()
    b.runTest()
