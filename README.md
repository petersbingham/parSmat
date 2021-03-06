# parsmat
Python package to parametrise the multi-channel S-matrix using a pade approximation according to the technique presented in "S.A. Rakityansky P.O.G. Ogunbade. S-matrix parametrization as a way of locating quantum resonances and bound states:multichannel case, 2010".

## Installation

Clone the repository and install with the following commands:

    git clone https://github.com/petersbingham/parsmat.git
    cd parsmat
    python setup.py install
    
## Dependencies
Author packages (these will have their own dependencies):
 - [pynumwrap](https://github.com/petersbingham/pynumwrap)
 - [channelutil](https://github.com/petersbingham/channelutil)
 - [tisutil](https://github.com/petersbingham/tisutil) (optional)
 - [stelempy](https://github.com/petersbingham/stelempy) (optional)

## Usage

There are four main functions `calculate_coefficients`, `get_elastic_Fin_fun`,  `get_elastic_Smat_fun` and `get_elastic_Qmat_fun`. Each of the functions take an `channelutil.AsymCalc` like object to describe the channels (see [`channelutil`](https://github.com/petersbingham/channelutil) and example below for more details).

There are two types that the ukrmolmatreader is compatible with, standard python types and mpmath types. Python types is the default. To change to mpmath types call the module function `use_mpmath_types()`.

#### Calculation of coefficients

Coefficients are calculated using the following function:

`calculate_coefficients(smatdata, asymcalc)`

The `smatdata`, is a dictionary like (eg python `dict` or a `tisutil.dSmat`) of S-matrices keyed by energy that will be used as input to the fit routine. It returns `alphas` and `betas`.

#### Calculation of scattering quantities

Once the coefficients have been calculated they can be passed to the following functions for calculation of scattering quantities. All these functions return either an energy function reference or the indicated `tisutil` like container.

`get_elastic_Fin_fun(coeffs, asymcalc)` returns `tisutil.cFinMatSympypolyk`.

`get_elastic_Smat_fun(coeffs, asymcalc)` returns `tisutil.cSmat`.

`get_elastic_Spmat_fun` returns `tisutil.cMat`

`get_elastic_Qmat_fun` returns `tisutil.cQmat`

The following example illustrates. Explanation follows.
```python
import channelutil as cu
import twochanradialwell as tcrw
import parsmat as psm

calc = cu.AsymCalc(cu.hartrees, [0,0])

cSmat = tcrw.get_Smat_fun(1.0,2.0,2.0,calc,1.0)
dSmat = cSmat.discretise(1.,8.,12)

coeffs = psm.calculate_coefficients(dSmat, calc)
cSmat = psm.get_elastic_Smat_fun(coeffs, calc)
dSmat = cSmat.discretise(1.,8.,100)
dSmat.plot()
```

The first code block imports the required dependent python packages. We are using [`twochanradialwell`](https://github.com/petersbingham/twochanradialwell) to provide our example fit data but the user can of course provide their own. The next single single code block prepares the asymptotic calculator, describing the channels of our calculation. We will look at an elastic `twochanradialwell`, so the second parameter specifies two channels, each with zero angular momentum. The next block of code prepares our demo fit data from the `twochanradialwell`, the second line here uses the `tisutil.dMat` interface to create the dictionary of S-matrices. In the final block of code the fit coefficients are calculated and then used to obtain a function reference to a parameterised S-matrix. Finally this fitted S-matrix is plotted using the `tisutil.dMat` interface.

The reference quoted at the beginning of this document describes a technique to find the poles of the S-matrix by looking for stable roots of the Fin. For this we recommend using the `find_roots` of the `cPolyVal` returned from `get_elastic_Fin_fun` (if tisutil is available) for several numbers of fit points and then using the `stelempy` package to locate the stable roots (corresponding to the S-matrix poles). Alternatively, the [`reskit`](https://github.com/petersbingham/reskit) package provides a tool combining all of this functionality.
