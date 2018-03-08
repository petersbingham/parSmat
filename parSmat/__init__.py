import pynumwrap as nw

import collections

class parSmatException(Exception):
    def __init__(self, string):
        self.string = string
    def __str__(self):
        return "Rad Well Error: " + self.string

def calculateCoefficients(sMatData, chanCalc):
    enes = [ene for ene in sorted(sMatData.keys(), key=lambda val: val.real)]
    _checkCoeffInput(enes, sMatData)
    return _calculateCoefficients(enes, sMatData, chanCalc)

def getFinFun(coeffs, chanCalc):
    pass

def getSmatFun(coeffs, chanCalc):
    pass


########################################################################   
####################### calculateCoefficients ##########################
########################################################################

def _checkCoeffInput(enes, sMatData):
    firstShape = nw.shape(sMatData[enes[0]])
    if firstShape[0]==0 or firstShape[0]!=firstShape[1]:
        raise parSmatException("Bad Input: Matrix not square")
    if firstShape != nw.shape(sMatData[enes[-1]]):
        raise parSmatException("S-matrices have difference shapes")
    if len(sMatData)<2:
        raise parSmatException("Bad Input: Not enough Data")
    if len(sMatData)<2:
        raise parSmatException("Bad Input: Specified fit size too small")
    if len(sMatData)>len(sMatData):
        raise parSmatException("Bad Input: Specified fitsize larger than number of data")
    if len(sMatData)%len(sMatData)!=0:
        raise parSmatException("Bad Input: Num of data is not a multiple of the fit size")
    if len(sMatData)%2!=0:
        raise parSmatException("Bad Input: Number of data not even")
    if len(sMatData)%2!=0:
        raise parSmatException("Bad Input: Fit size not even")

def _calculateCoefficients(enes, sMatData, chanCalc):
    numData = len(sMatData)
    numPolyTerms = numData / 2
    numCoeffs = numPolyTerms + 1
    numChannels = nw.shape(sMatData[enes[0]])[0]
    alphas = _initialiseCoefficients(numCoeffs, numChannels)
    betas = _initialiseCoefficients(numCoeffs, numChannels)

    for n in range(numChannels):
        sysMat = nw.matrix(_getSysMatInit(numData, numChannels, numPolyTerms))
        resVec = nw.matrix(_getResVecInit(numData, numChannels))
        for m in range(numChannels):
            ei = 0
            for ene in enes:
                for ti in range(numPolyTerms):  #We have two indices ci (coefficient) and ti (term). We know the first term in the poly expansion so numCoeffs = numPolyTerms + 1 
                    exp = ti+1
                    for j in range(numChannels): 
                        if j==m:
                            alphaCoeff = _primaryAlpha(sMatData, chanCalc, m, n, ene, exp)
                            betaCoeff = _primaryBeta(sMatData, chanCalc, m, n, ene, exp)
                        else:
                            alphaCoeff = _secondaryAlpha(sMatData, chanCalc, m, n, j, ene, exp)
                            betaCoeff = _secondaryBeta(sMatData, chanCalc, m, n, j, ene, exp)
                        sysMat[_row(numData,m,ei),_alphaIndex(numPolyTerms,j,ti)] = alphaCoeff
                        sysMat[_row(numData,m,ei),_betaIndex(numPolyTerms,numChannels,j,ti)] = betaCoeff
                resVec[_row(numData,m,ei),0] = _result(sMatData, m, n, ene)
                ei += 1
        coeffVec = nw.lin_solve(sysMat, resVec)
        _copyColumnCoeffs(alphas, betas, coeffVec, numPolyTerms, numChannels, numCoeffs, n)
    return alphas, betas

def _initialiseCoefficients(numCoeffs, numChannels):
    coeffs = []
    for _ in range(0, numCoeffs):
        mat = nw.matrix(_getZeroListMats(numChannels))
        coeffs.append(mat)
    return coeffs

def _getZeroListMats(numChannels):
    return [[0.0+0.0j]*numChannels]*numChannels

def _getSysMatInit(numData, numChannels, numPolyTerms):
    return [[0.0]*2*numPolyTerms*numChannels]*numData*numChannels
   
def _getResVecInit(numData, numChannels):
    return [[0.0]]*numData*numChannels


def _primaryAlpha(sMatData, chanCalc, m, n, ene, exp):
    return _kl(chanCalc,n,ene,1.0) / _kl(chanCalc,m,ene,1.0) * (sMatData[ene][m,m]-1.0) * nw.pow(ene,exp)

def _primaryBeta(sMatData, chanCalc, m, n, ene, exp):
    return -1.0j * _kl(chanCalc,m,ene,0.0) * _kl(chanCalc,n,ene,1.0) * (sMatData[ene][m,m]+1.0) * nw.pow(ene,exp)

def _secondaryAlpha(sMatData, chanCalc, m, n, j, ene, exp):
    return _kl(chanCalc,n,ene,1.0) / _kl(chanCalc,j,ene,1.0) * sMatData[ene][m,j] * nw.pow(ene,exp)

def _secondaryBeta(sMatData, chanCalc, m, n, j, ene, exp):
    return -1.0j * _kl(chanCalc,j,ene,0.0) * _kl(chanCalc,n,ene,1.0) * sMatData[ene][m,j] * nw.pow(ene,exp)


def _row(numData, m, ei):
    return m*numData + ei

def _alphaIndex(numPolyTerms, m, ti):
    return m*numPolyTerms + ti

def _betaIndex(numPolyTerms, numChannels, m, ti):
    return numPolyTerms*numChannels + m*numPolyTerms + ti


def _result(sMatData, m, n, ene):
    num = 0.0
    if m==n:
        num = 1.0
    return num - sMatData[ene][m,n]

def _copyColumnCoeffs(alphas, betas, coeffVec, numPolyTerms, numChannels, numCoeffs, n):
    for ci in range(numCoeffs):
        ti = ci-1
        for m in range(numChannels):
            if ci==0:
                if m==n:
                    alphas[ci][m,n] = 1.0
            else:
                alphas[ci][m,n] = nw.complex(coeffVec[_alphaIndex(numPolyTerms,m,ti),0])
                betas[ci][m,n] = nw.complex(coeffVec[_betaIndex(numPolyTerms,numChannels,m,ti),0])

def _kl(chanCalc, ch, ene, mod):
    k = chanCalc.k(ch, ene)
    return nw.pow(k, chanCalc.ls[ch]+mod)

def usePythonTypes(dps=nw.dps_default_python):
    nw.usePythonTypes(dps)

def useMpmathTypes(dps=nw.dps_default_mpmath):
    nw.useMpmathTypes(dps)
