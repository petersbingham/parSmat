import pynumwrap as nw
import sympy as sym
from sympy.matrices import Matrix as sym_matrix
import collections
try:
    import tisutil as tu
except:
    tu = None

########################################################################   
###################### Calculate Coefficients ##########################
########################################################################

def _checkCoeffInput(enes, sMatData, asymCal):
    firstShape = nw.shape(sMatData[enes[0]])
    excStr = ""
    if firstShape[0]==0 or firstShape[0]!=firstShape[1]:
        excStr = "Bad Input: Matrix not square"
    if firstShape[0]!=asymCal.getNumberChannels():
        excStr = "Bad Input: Inconsistent channel specification"
    if firstShape != nw.shape(sMatData[enes[-1]]):
        excStr = "Bad Input: S-matrices have difference shapes"
    if len(sMatData)<2:
        excStr = "Bad Input: Not enough Data"
    if len(sMatData)<2:
        excStr = "Bad Input: Specified fit size too small"
    if len(sMatData)>len(sMatData):
        excStr = "Bad Input: Specified fitsize larger than number of data"
    if len(sMatData)%len(sMatData)!=0:
        excStr = "Bad Input: Num of data is not a multiple of the fit size"
    if len(sMatData)%2!=0:
        excStr = "Bad Input: Number of data not even"
    if len(sMatData)%2!=0:
        excStr = "Bad Input: Fit size not even"
    if excStr != "":
        raise parSmatException(excStr)

def _calculateCoefficients(enes, sMatData, asymCal):
    numData = len(sMatData)
    if not asymCal.isPolar():
        numUnknownTerms = numData / 2
        numCoeffs = numUnknownTerms + 1
    else:
        numUnknownTerms = numData / 2
        numCoeffs = numUnknownTerms
    numChannels = nw.shape(sMatData[enes[0]])[0]
    alphas = _initialiseCoefficients(numCoeffs, numChannels)
    betas = _initialiseCoefficients(numCoeffs, numChannels)

    for n in range(numChannels):
        sysMat = nw.matrix(_getSysMatInit(numData, numChannels, 
                                          numUnknownTerms))
        if not asymCal.isPolar():
            resVec = nw.matrix(_getResVecInit(numData, numChannels))
        for m in range(numChannels):
            ei = 0
            for ene in enes:
                # Two indices ci (coefficient) and ti (term).
                for ti in range(numUnknownTerms):
                    if not asymCal.isPolar():
                        # We know the first term in the poly expansion.
                        exp = ti+1
                    else:
                        exp = ti
                    for j in range(numChannels): 
                        if j==m:
                            alphaCoeff = _primaryAlpha(sMatData, asymCal, 
                                                       m, n, ene, exp)
                            betaCoeff = _primaryBeta(sMatData, asymCal, 
                                                     m, n, ene, exp)
                        else:
                            alphaCoeff = _secondaryAlpha(sMatData, asymCal, 
                                                         m, n, j, ene, exp)
                            betaCoeff = _secondaryBeta(sMatData, asymCal, 
                                                       m, n, j, ene, exp)
                        r = _row(numData,m,ei)
                        c = _alphaIndex(numUnknownTerms,j,ti)
                        sysMat[r,c] = alphaCoeff
                        c = _betaIndex(numUnknownTerms, numChannels, j, ti)
                        sysMat[r,c] = betaCoeff
                if not asymCal.isPolar():
                    resVec[_row(numData,m,ei),0] = _result(sMatData, m, n, ene)
                ei += 1
        if not asymCal.isPolar():
            coeffVec = nw.lin_solve(sysMat, resVec)
        else:
            coeffVec = nw.lin_solve_homo(sysMat)
        _copyColumnCoeffs(alphas, betas, coeffVec, asymCal, numUnknownTerms, 
                          numChannels, numCoeffs, n)
    return alphas, betas

def _initialiseCoefficients(numCoeffs, numChannels):
    coeffs = []
    for _ in range(0, numCoeffs):
        mat = nw.matrix(_getZeroListMats(numChannels))
        coeffs.append(mat)
    return coeffs

def _getZeroListMats(numChannels):
    return [[0.0+0.0j]*numChannels]*numChannels

def _getSysMatInit(numData, numChannels, numUnknownTerms):
    return [[0.0]*2*numUnknownTerms*numChannels]*numData*numChannels
   
def _getResVecInit(numData, numChannels):
    return [[0.0]]*numData*numChannels


def _primaryAlpha(sMatData, asymCal, m, n, ene, exp):
    kcal = _kl(asymCal,n,ene,1.0) / _kl(asymCal,m,ene,1.0)
    return kcal * (sMatData[ene][m,m]-1.0) * nw.pow(ene,exp)

def _primaryBeta(sMatData, asymCal, m, n, ene, exp):
    kcal = _kl(asymCal,m,ene,0.0) * _kl(asymCal,n,ene,1.0)
    return -1.0j * kcal * (sMatData[ene][m,m]+1.0) * nw.pow(ene,exp)

def _secondaryAlpha(sMatData, asymCal, m, n, j, ene, exp):
    kcal = _kl(asymCal,n,ene,1.0) / _kl(asymCal,j,ene,1.0)
    return kcal * sMatData[ene][m,j] * nw.pow(ene,exp)

def _secondaryBeta(sMatData, asymCal, m, n, j, ene, exp):
    kcal = _kl(asymCal,j,ene,0.0) * _kl(asymCal,n,ene,1.0)
    return -1.0j * kcal * sMatData[ene][m,j] * nw.pow(ene,exp)


def _row(numData, m, ei):
    return m*numData + ei

def _alphaIndex(numUnknownTerms, m, ti):
    return m*numUnknownTerms + ti

def _betaIndex(numUnknownTerms, numChannels, m, ti):
    return numUnknownTerms*numChannels + m*numUnknownTerms + ti


def _result(sMatData, m, n, ene):
    num = 0.0
    if m==n:
        num = 1.0
    return num - sMatData[ene][m,n]

def _copyColumnCoeffs(alphas, betas, coeffVec, asymCal, numUnknownTerms,
                      numChannels, numCoeffs, n):
    for ci in range(numCoeffs):
        if not asymCal.isPolar():
            ti = ci-1
        else:
            ti = ci
        for m in range(numChannels):
            if not asymCal.isPolar() and ci==0:
                if m==n:
                    alphas[ci][m,n] = 1.0
            else:
                r = _alphaIndex(numUnknownTerms,m,ti)
                alphas[ci][m,n] = nw.complex(coeffVec[r,0])
                r = _betaIndex(numUnknownTerms,numChannels,m,ti)
                betas[ci][m,n] = nw.complex(coeffVec[r,0])

def _kl(asymCal, ch, ene, mod):
    k = asymCal.k(ch, ene)
    return nw.pow(k, asymCal.ls[ch]+mod)

########################################################################   
###################### Parameterised Functions #########################
########################################################################

def _getElasticMatrix(coeffs, asymCal, finOnly, **kwargs):
    alphas = coeffs[0]
    betas = coeffs[1]
    numChannels = asymCal.getNumberChannels()
    k = nw.sym.symbols('k')
    matLst_fin = []
    if not finOnly:
        matLst_fout = []
    fact1 = (1.0/2.0)
    fact2 = 1.0/asymCal.getEneConv()
    for m in range(numChannels):
        matLst_fin.append([])
        if not finOnly:
            matLst_fout.append([])
        for n in range(numChannels):
            lm = asymCal.l(m)
            ln = asymCal.l(n)
            val_fin = 0.0
            val_fout = 0.0
            for ci in range(len(coeffs[0])):
                A = alphas[ci][m,n]
                B = betas[ci][m,n]
                real = nw.toSympy(A)*k**(ln-lm+2*ci)
                imag = sym.I*nw.toSympy(B)*k**(ln+lm+1+2*ci)
                fact3 = fact1*fact2**ci
                v = fact3 * (real - imag)
                val_fin += v
                if not finOnly:
                    val_fout += fact3 * (real + imag)
            matLst_fin[len(matLst_fin)-1].append(val_fin)
            if not finOnly:
                matLst_fout[len(matLst_fout)-1].append(val_fout)
    mat_fin = sym_matrix(matLst_fin)
    if not finOnly:
        mat_fout = sym_matrix(matLst_fout)
        if "sym_matrix_inv" in kwargs:
            mat_fin_inv = mat_fin.inv(**kwargs["sym_matrix_inv"])
        else:
            mat_fin_inv = mat_fin.inv()
        return mat_fout * mat_fin_inv
    else:
        return mat_fin

########################################################################   
######################### Public Interface #############################
########################################################################

class parSmatException(Exception):
    def __init__(self, string):
        self.string = string
    def __str__(self):
        return "Rad Well Error: " + self.string

def calculateCoefficients(sMatData, asymCal):
    enes = [ene for ene in sorted(sMatData.keys(), key=lambda val: val.real)]
    _checkCoeffInput(enes, sMatData, asymCal)
    return _calculateCoefficients(enes, sMatData, asymCal)

def getElasticFinFun(coeffs, asymCal):
    mat = _getElasticMatrix(coeffs, asymCal, True)
    ret = lambda ene: nw.fromSympyMatrix(mat.subs('k', asymCal.fk(ene)))
    if tu is not None:
        ret = tu.cPolykmat(mat, 'k', asymCal)
    return ret

def getElasticSmatFun(coeffs, asymCal, **kwargs):
    mat = _getElasticMatrix(coeffs, asymCal, False, **kwargs)
    ret = lambda ene: nw.fromSympyMatrix(mat.subs('k', asymCal.fk(ene)))
    if tu is not None:
        ret = tu.cPolySmat(mat, 'k', asymCal)
    return ret

# Ancillary helper functions:
def getNumCoeffForN(N):
    return N/2 + 1

# Type functions:
def usePythonTypes(dps=nw.dps_default_python):
    nw.usePythonTypes(dps)

def useMpmathTypes(dps=nw.dps_default_mpmath):
    nw.useMpmathTypes(dps)

def setTypeMode(mode, dps=None):
    nw.setTypeMode(mode, dps)
