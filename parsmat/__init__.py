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
    first_shape = nw.shape(smatdata[enes[0]])
    excStr = ""
    if first_shape[0]==0 or first_shape[0]!=first_shape[1]:
        excStr = "Bad Input: Matrix not square"
    if first_shape[0]!=asymcalc.get_number_channels():
        excStr = "Bad Input: Inconsistent channel specification"
    if first_shape != nw.shape(smatdata[enes[-1]]):
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
    num_data = len(smatdata)
    num_poly_terms = num_data / 2
    num_coeffs = num_poly_terms + 1
    num_channels = nw.shape(smatdata[enes[0]])[0]
    alphas = _initialise_coefficients(num_coeffs, num_channels)
    betas = _initialise_coefficients(num_coeffs, num_channels)

    for j in range(num_channels):
        sys_mat = nw.matrix(_get_sys_mat_init(num_data, num_channels, num_poly_terms))
        res_vec = nw.matrix(_get_res_vec_init(num_data, num_channels))
        for i in range(num_channels):
            ei = 0
            for ene in enes:
                for ti in range(num_poly_terms):  #We have two indices ci (coefficient) and ti (term). We know the first term in the poly expansion so num_coeffs = num_poly_terms + 1 
                    exp = ti+1
                    for k in range(num_channels): 
                        if k==i:
                            alpha_coeff = _primary_alpha(smatdata, asymcalc, 
                                                         i, j, ene, exp)
                            beta_coeff = _primary_beta(smatdata, asymcalc, 
                                                       i, j, ene, exp)
                        else:
                            alpha_coeff = _secondary_alpha(smatdata, asymcalc, 
                                                           i, j, k, ene, exp)
                            beta_coeff = _secondary_beta(smatdata, asymcalc, 
                                                         i, j, k, ene, exp)
                        sys_mat[_row(num_data,i,ei),_alpha_index(num_poly_terms,k,ti)] = alpha_coeff
                        sys_mat[_row(num_data,i,ei),_beta_index(num_poly_terms,num_channels,k,ti)] = beta_coeff
                res_vec[_row(num_data,i,ei),0] = _result(smatdata, i, j, ene)
                ei += 1
        coeff_vec = nw.lin_solve(sys_mat, res_vec)
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

def _get_sys_mat_init(num_data, num_channels, num_poly_terms):
    return [[0.0]*2*num_poly_terms*num_channels]*num_data*num_channels
   
def _get_res_vec_init(num_data, num_channels):
    return [[0.0]]*num_data*num_channels


def _primary_alpha(smatdata, asymcalc, i, j, ene, exp):
    return _kl(asymcalc,j,ene,1.0) / _kl(asymcalc,i,ene,1.0) * (smatdata[ene][i,i]-1.0) * nw.pow(ene,exp)

def _primary_beta(smatdata, asymcalc, i, j, ene, exp):
    return -1.0j * _kl(asymcalc,i,ene,0.0) * _kl(asymcalc,j,ene,1.0) * (smatdata[ene][i,i]+1.0) * nw.pow(ene,exp)

def _secondary_alpha(smatdata, asymcalc, i, j, k, ene, exp):
    return _kl(asymcalc,j,ene,1.0) / _kl(asymcalc,k,ene,1.0) * smatdata[ene][i,k] * nw.pow(ene,exp)

def _secondary_beta(smatdata, asymcalc, i, j, k, ene, exp):
    return -1.0j * _kl(asymcalc,k,ene,0.0) * _kl(asymcalc,j,ene,1.0) * smatdata[ene][i,k] * nw.pow(ene,exp)


def _row(num_data, i, ei):
    return i*num_data + ei

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
    return nw.pow(k, asymcalc.angmoms[ch]+mod)

########################################################################   
###################### Parameterised Functions #########################
########################################################################

def _get_elastic_matrix(coeffs, asymcalc, calc, k):
    alphas = coeffs[0]
    betas = coeffs[1]
    num_channels = asymcalc.get_number_channels()
    fact1 = 1.0/2.0
    fact2 = 1.0/asymcalc.get_ene_conv()
    for i in range(num_channels):
        calc.chan_start()
        for j in range(num_channels):
            lm = asymcalc.angmom(i)
            ln = asymcalc.angmom(j)
            calc.reset_terms()
            for ci in range(len(coeffs[0])):
                A = alphas[ci][i,j]
                B = betas[ci][i,j]
                calc.calc_term(A, B, k, ln, lm, ci, fact1*fact2**ci)
            calc.set_terms()
    calc.finialise()
    return calc.get_matrix()

## Matrix Calculators

class Fin_calc:
    def __init__(self):
        self.list_fin = []
        self.mat_fin = None
    def chan_start(self):
        self.list_fin.append([])
    def reset_terms(self):
        self.val_fin = 0.
    def calc_term(self, A, B, k, ln, lm, ci, fact):
        real = nw.to_sympy(A)*k**(ln-lm+2*ci)
        imag = nw.to_sympy(B)*sym.I*k**(ln+lm+1+2*ci)
        self.val_fin += fact * (real - imag)
    def set_terms(self):
        self.list_fin[len(self.list_fin)-1].append(self.val_fin)
    def finialise(self):
        self.mat_fin = sym_matrix(self.list_fin)
    def get_matrix(self):
        return self.mat_fin

class S_calc:
    def __init__(self):
        self.list_fin = []
        self.list_fout = []
        self.mat_fin = None
        self.mat_fout = None
    def chan_start(self):
        self.list_fin.append([])
        self.list_fout.append([])
    def reset_terms(self):
        self.val_fin = 0.
        self.val_fout = 0.
    def calc_term(self, A, B, k, ln, lm, ci, fact):
        real = A*k**(ln-lm+2*ci)
        imag = B*1.j*k**(ln+lm+1+2*ci)
        self.val_fin += fact * (real - imag)
        self.val_fout += fact * (real + imag)
    def set_terms(self):
        self.list_fin[len(self.list_fin)-1].append(self.val_fin)
        self.list_fout[len(self.list_fout)-1].append(self.val_fout)
    def finialise(self):
        self.mat_fin = nw.matrix(self.list_fin)
        self.mat_fout = nw.matrix(self.list_fout)
    def get_matrix(self):
        return self.mat_fout * nw.invert(self.mat_fin)

class Sp_calc:
    def __init__(self, s_calc):
        self.s_calc = s_calc
        self.list_fin_p = []
        self.list_fout_p = []
        self.mat_fin_p = []
        self.mat_fout_p = []
    def chan_start(self):
        self.s_calc.chan_start()
        self.list_fin_p.append([])
        self.list_fout_p.append([])
    def reset_terms(self):
        self.s_calc.reset_terms()
        self.val_fin_p = 0.
        self.val_fout_p = 0.
    def calc_term(self, A, B, k, ln, lm, ci, fact):
        self.s_calc.calc_term(A, B, k, ln, lm, ci, fact)
        real = A*k**(ln-lm+2*ci-2) * (ln-lm+2*ci)
        imag = B*1.j*k**(ln+lm-1+2*ci) * (ln+lm+1+2*ci)
        self.val_fin_p += (1./2.) * fact * (real - imag)
        self.val_fout_p += (1./2.) * fact * (real + imag)
    def set_terms(self):
        self.s_calc.set_terms()
        self.list_fin_p[len(self.list_fin_p)-1].append(self.val_fin_p)
        self.list_fout_p[len(self.list_fout_p)-1].append(self.val_fout_p)
    def finialise(self):
        self.s_calc.finialise()
        self.mat_fin_p = nw.matrix(self.list_fin_p)
        self.mat_fout_p = nw.matrix(self.list_fout_p)
    def get_matrix(self):
        fin = self.s_calc.mat_fin
        fin_inv = nw.invert(fin)
        fin_p = self.mat_fin_p
        fin_inv_p = -fin_inv * fin_p * fin_inv

        fout = self.s_calc.mat_fout
        fout_p = self.mat_fout_p

        return fout * fin_inv_p  +  fout_p * fin_inv

class Q_calc:
    def __init__(self):
        self.s_calc = S_calc()
        self.sp_calc = Sp_calc(self.s_calc)
    def chan_start(self):
        self.sp_calc.chan_start()
    def reset_terms(self):
        self.sp_calc.reset_terms()
    def calc_term(self, A, B, k, ln, lm, ci, fact):
        self.sp_calc.calc_term(A, B, k, ln, lm, ci, fact)
    def set_terms(self):
        self.sp_calc.set_terms()
    def finialise(self):
        self.sp_calc.finialise()
    def get_matrix(self):
        smat = self.s_calc.get_matrix()
        smat_p = self.sp_calc.get_matrix()
        return 1.j*smat*nw.dagger(smat_p)

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
    mat = _get_elastic_matrix(coeffs, asymcalc, Fin_calc(), nw.sym.symbols('k'))
    ret = lambda ene: nw.from_sympy_matrix(mat.subs('k', asymcalc.fk(ene)))
    if tu is not None:
        ret = tu.cFinMatSympypolyk(mat, 'k', asymcalc)
    return ret

def get_elastic_Smat_fun(coeffs, asymcalc):
    ret = lambda ene: _get_elastic_matrix(coeffs, asymcalc, S_calc(),
                                          asymcalc.fk(ene))
    if tu is not None:
        ret = tu.cSmat(ret, asymcalc)
    return ret

def get_elastic_Spmat_fun(coeffs, asymcalc):
    ret = lambda ene: _get_elastic_matrix(coeffs, asymcalc, Sp_calc(S_calc()),
                                          asymcalc.fk(ene))
    if tu is not None:
        ret = tu.cMat(ret, asymcalc)
    return ret

def get_elastic_Qmat_fun(coeffs, asymcalc):
    ret = lambda ene: _get_elastic_matrix(coeffs, asymcalc, Q_calc(),
                                          asymcalc.fk(ene))
    if tu is not None:
        ret = tu.cQmat(ret, asymcalc)
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
