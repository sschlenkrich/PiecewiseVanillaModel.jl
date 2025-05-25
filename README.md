# PiecewiseVanillaModel.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sschlenkrich.github.io/PiecewiseVanillaModel.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sschlenkrich.github.io/PiecewiseVanillaModel.jl/dev/)
[![Build Status](https://github.com/sschlenkrich/PiecewiseVanillaModel.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sschlenkrich/PiecewiseVanillaModel.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/sschlenkrich/PiecewiseVanillaModel.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sschlenkrich/PiecewiseVanillaModel.jl)


The PiecewiseVanillaModel.jl package implements a static model for the pricing and implied volatility calculation of Vanilla options.

The model is based on a piece-wise specification of the terminal distribution of the underlying financial risk factor.

Details of the model are discussed in [S. Schlenkrich, A Piece-wise Model for Vanilla Option Pricing](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5214566).

## Installation

The package can be installed from this Github repository.

```
using Pkg
Pkg.add(url="https://github.com/sschlenkrich/PiecewiseVanillaModel.jl")
```

## Usage

Suppose we have a given set of implied Black volatilities for a risk factor with forward price $S_0 = 1.00$ and time to option exercise $T=2.0$.

```
S0 = 1.00
T = 2.0
absolute_strikes = [ 0.50, 0.75, 0.90, 1.10, 1.50, 2.00 ]

sigma_b76_atm = 0.30
sigma_b76_smile = [ 0.30, 0.30, 0.30, 0.30, 0.30, 0.30 ]
```

To calibrate a model, we need to convert Black volatilities to Normal volatilities. We opt for Normal volatilities as calibration inputs because they allow for negative risk factor values and strikes while being independent of option notional and call/put flag.

```
using PiecewiseVanillaModel
pvm = PiecewiseVanillaModel

sigma_n_atm = pvm.normal_volatility(sigma_b76_atm, S0, S0, T)
sigma_n_smile = [
    pvm.normal_volatility(v, K, S0, T)
    for (v, K) in zip(sigma_b76_smile, absolute_strikes)
]
```

Moreover, we want to express strikes as relative strikes.

```
relative_strikes = absolute_strikes .- S0
```

A key modelling input is the specification of the reference strikes. The reference strikes are specified in terms of *standard deviations from ATM*; assuming constant normal volatility.

```
stdevs_lower_smile = [ 0.5, 1.0 ]
stdevs_upper_smile = [ 0.5, 1.0 ]
```

Furthermore, we specify smile extrapolation using linear extrapolation of volatility parameters.

```
r_extrap_lower_smile = nothing
r_extrap_upper_smile = nothing
```

Finally, we can construct and calibrate a model.

```
res = pvm.calibrated_model_from_smile(
    S0, sigma_n_atm, T,
    stdevs_lower_smile, stdevs_upper_smile,
    relative_strikes, sigma_n_smile,
    rexl = r_extrap_lower_smile,
    rexu = r_extrap_upper_smile,
)
```

The function returns a named tuple with attributes `model` and `result`. `model` represents the calibrated model and `result` contains calibration result details of type `LsqFit.LsqFitResult`.

We can calculate, for example, at-the-money call and put prices.

```
pvm.call_option(res.model, S0)
pvm.put_option(res.model, S0)
```

This should yield `0.167995971427363`.

Similarly, we can double-check the implied volatility for the forward (or any other strike).

```
pvm.lognormal_volatility(res.model, S0)
```

This returns `0.30`, which was also an input to the calibration.

It turns out that the calibrated model in this example is a log-normal model with log-normal volatility of 30%. This can be verified by inspecting the model attributes.

```
res.model.v0
# 0.3000000000057258

res.model.w0
# 0.29999999999865457

res.model.dvu ./ res.model.dsu
# 2-element Vector{Float64}:
#  0.29999999994989296
#  0.29999999997949833

-res.model.dvl ./ res.model.dsl
# 2-element Vector{Float64}:
#  0.3000000000402363
#  0.300000000014485
```

Model usage is also illustrated in the example Julia notebook [JaeckelExample.ipynb](./notebooks/JaeckelExample.ipynb).

An alternative Python implementation with an example is provided in [PiecewiseVanillaModel.ipynb](./python/PiecewiseVanillaModel.ipynb).

