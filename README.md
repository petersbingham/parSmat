# parsmat
Python package to parametrise the multi-channel S-matrix using a pade approximation according to the technique presented in "S.A. Rakityansky P.O.G. Ogunbade. S-matrix parametrization as a way of locating quantum resonances and bound states:multichannel case, 2010".

## Installation

Clone the repository and install with the following commands:

    git clone https://github.com/petersbingham/parsmat.git
    cd parsmat
    python setup.py install
    
## Dependencies
Author Libraries (these will have their own dependencies):
 - pynumwrap https://github.com/petersbingham/pynumwrap
 - channelutil (recommended) https://github.com/petersbingham/channelutil
 - tisutil (optional) https://github.com/petersbingham/tisutil

## Usage

There are three main functions `calculate_coefficients`, `get_elastic_Fin_fun` and `get_elastic_Smat_fun`. Each of the functions take an `channelutil.AsymCalc` like object to describe the channels (see https://github.com/petersbingham/channelutil and example below for more details).

There are two types that the ukrmolmatreader is compatible with, standard python types and mpmath types. Python types is the default. To change to mpmath types call the module function `use_mpmath_types()`.

#### `calculate_coefficients(smatdata, asymcalc)`

The `smatdata`, is a dictionary like (eg python `dict` or a `tisutil.dSmat`) of S-matrices keyed by energy that will be used as input to the fit routine. It returns `alphas` and `betas`.

#### `get_elastic_Fin_fun(coeffs, asymcalc)`

Using the `alphas` and `betas` coefficients returned from `calculate_coefficients` this function returns either an energy function reference or a `tisutil.cPolykmat` describing the parameterised Jost denominator (Fin).

#### `get_elastic_Smat_fun(coeffs, asymcalc)`

Using the `alphas` and `betas` coefficients returned from `calculate_coefficients` this function returns either an energy function reference or a `tisutil.cSmat` describing the parameterised S-matrix.

The following example illustrates. Explanation follows.
```python
>>> import parsmat as pm
```
