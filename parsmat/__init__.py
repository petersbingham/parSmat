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

def _check_coeff_input(enes, smatdata, asymcalc):
    firstShape = nw.shape(smatdata[enes[0]])
    excStr = ""
    if firstShape[0]==0 or firstShape[0]!=firstShape[1]:
        excStr = "Bad Input: Matrix not square"
    if firstShape[0]!=asymcalc.get_number_channels():
        excStr = "Bad Input: Inconsistent channel specification"
    if firstShape != nw.shape(smatdata[enes[-1]]):
        excStr = "Bad Input: S-matrices have difference shapes"
    if len(smatdata)<2:
        excStr = "Bad Input: Not enough Data"
    if len(smatdata)<2:
        excStr = "Bad Input: Specified fit size too small"
    if len(smatdata)>len(smatdata):
        excStr = "Bad Input: Specified fitsize larger than number of data"
    if len(smatdata)%len(smatdata)!=0:
        excStr = "Bad Input: Num of data is not a multiple of the fit size"
    if len(smatdata)%2!=0:
        excStr = "Bad Input: Number of data not even"
    if len(smatdata)%2!=0:
        excStr = "Bad Input: Fit size not even"
    if excStr != "":
        raise ParSmatException(excStr)

def _calculate_coefficients(enes, smatdata, asymcalc):
    numData = len(smatdata)
    num_poly_terms = numData / 2
    num_coeffs = num_poly_terms + 1
    num_channels = nw.shape(smatdata[enes[0]])[0]
    alphas = _initialise_coefficients(num_coeffs, num_channels)
    betas = _initialise_coefficients(num_coeffs, num_channels)

    for j in range(num_channels):
        sysMat = nw.matrix(_get_sys_mat_init(numData, num_channels, num_poly_terms))
        resVec = nw.matrix(_get_res_vec_init(numData, num_channels))
        for i in range(num_channels):
            ei = 0
            for ene in enes:
                for ti in range(num_poly_terms):  #We have two indices ci (coefficient) and ti (term). We know the first term in the poly expansion so num_coeffs = num_poly_terms + 1 
                    exp = ti+1
                    for k in range(num_channels): 
                        if k==i:
                            alphaCoeff = _primary_alpha(smatdata, asymcalc, 
                                                       i, j, ene, exp)
                            betaCoeff = _primary_beta(smatdata, asymcalc, 
                                                     i, j, ene, exp)
                        else:
                            alphaCoeff = _secondary_alpha(smatdata, asymcalc, 
                                                         i, j, k, ene, exp)
                            betaCoeff = _secondary_beta(smatdata, asymcalc, 
                                                       i, j, k, ene, exp)
                        sysMat[_row(numData,i,ei),_alpha_index(num_poly_terms,k,ti)] = alphaCoeff
                        sysMat[_row(numData,i,ei),_beta_index(num_poly_terms,num_channels,k,ti)] = betaCoeff
                resVec[_row(numData,i,ei),0] = _result(smatdata, i, j, ene)
                ei += 1
        coeff_vec = nw.lin_solve(sysMat, resVec)
        _copy_column_coeffs(alphas, betas, coeff_vec, num_poly_terms, num_channels, num_coeffs, j)
    return alphas, betas

def _initialise_coefficients(num_coeffs, num_channels):
    coeffs = []
    for _ in range(0, num_coeffs):
        mat = nw.matrix(_get_zero_list_mats(num_channels))
        coeffs.append(mat)
    return coeffs

def _get_zero_list_mats(num_channels):
    return [[0.0+0.0j]*num_channels]*num_channels

def _get_sys_mat_init(numData, num_channels, num_poly_terms):
    return [[0.0]*2*num_poly_terms*num_channels]*numData*num_channels
   
def _get_res_vec_init(numData, num_channels):
    return [[0.0]]*numData*num_channels


def _primary_alpha(smatdata, asymcalc, i, j, ene, exp):
    return _kl(asymcalc,j,ene,1.0) / _kl(asymcalc,i,ene,1.0) * (smatdata[ene][i,i]-1.0) * nw.pow(ene,exp)

def _primary_beta(smatdata, asymcalc, i, j, ene, exp):
    return -1.0j * _kl(asymcalc,i,ene,0.0) * _kl(asymcalc,j,ene,1.0) * (smatdata[ene][i,i]+1.0) * nw.pow(ene,exp)

def _secondary_alpha(smatdata, asymcalc, i, j, k, ene, exp):
    return _kl(asymcalc,j,ene,1.0) / _kl(asymcalc,k,ene,1.0) * smatdata[ene][i,k] * nw.pow(ene,exp)

def _secondary_beta(smatdata, asymcalc, i, j, k, ene, exp):
    return -1.0j * _kl(asymcalc,k,ene,0.0) * _kl(asymcalc,j,ene,1.0) * smatdata[ene][i,k] * nw.pow(ene,exp)


def _row(numData, i, ei):
    return i*numData + ei

def _alpha_index(num_poly_terms, i, ti):
    return i*num_poly_terms + ti

def _beta_index(num_poly_terms, num_channels, i, ti):
    return num_poly_terms*num_channels + i*num_poly_terms + ti


def _result(smatdata, i, j, ene):
    num = 0.0
    if i==j:
        num = 1.0
    return num - smatdata[ene][i,j]

def _copy_column_coeffs(alphas, betas, coeff_vec, num_poly_terms, num_channels, num_coeffs, j):
    for ci in range(num_coeffs):
        ti = ci-1
        for i in range(num_channels):
            if ci==0:
                if i==j:
                    alphas[ci][i,j] = 1.0
            else:
                alphas[ci][i,j] = nw.complex(coeff_vec[_alpha_index(num_poly_terms,i,ti),0])
                betas[ci][i,j] = nw.complex(coeff_vec[_beta_index(num_poly_terms,num_channels,i,ti),0])

def _kl(asymcalc, ch, ene, mod):
    k = asymcalc.k(ch, ene)
    return nw.pow(k, asymcalc.ls[ch]+mod)

########################################################################   
###################### Parameterised Functions #########################
########################################################################

def _convert(fin_only, val, imag=False):
    if fin_only:
        v = nw.to_sympy(val)
        if imag:
            v *= sym.I
    else:
        v = val
        if imag:
            v *= 1.j
    return v

def _get_elastic_matrix(coeffs, asymcalc, fin_only, k):
    alphas = coeffs[0]
    betas = coeffs[1]
    num_channels = asymcalc.get_number_channels()
    matLst_fin = []
    if not fin_only:
        matLst_fout = []
    fact1 = (1.0/2.0)
    fact2 = 1.0/asymcalc.get_ene_conv()
    for i in range(num_channels):
        matLst_fin.append([])
        if not fin_only:
            matLst_fout.append([])
        for j in range(num_channels):
            lm = asymcalc.l(i)
            ln = asymcalc.l(j)
            val_fin = 0.0
            val_fout = 0.0
            for ci in range(len(coeffs[0])):
                A = alphas[ci][i,j]
                B = betas[ci][i,j]
                real = _convert(fin_only,A)*k**(ln-lm+2*ci)
                imag = _convert(fin_only,B,True)*k**(ln+lm+1+2*ci)
                fact3 = fact1*fact2**ci
                v = fact3 * (real - imag)
                val_fin += v
                if not fin_only:
                    val_fout += fact3 * (real + imag)
            matLst_fin[len(matLst_fin)-1].append(val_fin)
            if not fin_only:
                matLst_fout[len(matLst_fout)-1].append(val_fout)
    if not fin_only:
        mat_fin = nw.matrix(matLst_fin)
        mat_fout = nw.matrix(matLst_fout)
        return mat_fout * nw.invert(mat_fin)
    else:
        return sym_matrix(matLst_fin)

########################################################################   
######################### Public Interface #############################
########################################################################

class ParSmatException(Exception):
    def __init__(self, string):
        self.string = string
    def __str__(self):
        return "Rad Well Error: " + self.string

def calculate_coefficients(smatdata, asymcalc):
    enes = [ene for ene in sorted(smatdata.keys(), key=lambda val: val.real)]
    _check_coeff_input(enes, smatdata, asymcalc)
    return _calculate_coefficients(enes, smatdata, asymcalc)

def get_elastic_Fin_fun(coeffs, asymcalc):
    mat = _get_elastic_matrix(coeffs, asymcalc, True, nw.sym.symbols('k'))
    ret = lambda ene: nw.from_sympy_matrix(mat.subs('k', asymcalc.fk(ene)))
    if tu is not None:
        ret = tu.cPolykmat(mat, 'k', asymcalc)
    return ret

def get_elastic_Smat_fun(coeffs, asymcalc):
    funref = lambda ene: _get_elastic_matrix(coeffs, asymcalc, False, 
                                           asymcalc.fk(ene))
    if tu is not None:
        ret = tu.cSmat(funref, asymcalc)
    return ret

# Ancillary helper functions:
def get_num_coeff_for_Npts(Npts):
    return Npts/2 + 1

# Type functions:
def use_python_types(dps=nw.dps_default_python):
    nw.use_python_types(dps)

def use_mpmath_types(dps=nw.dps_default_mpmath):
    nw.use_mpmath_types(dps)

def set_type_mode(mode, dps=None):
    nw.set_type_mode(mode, dps)
