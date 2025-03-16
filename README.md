# PiecewiseVanillaModel.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sschlenkrich.github.io/PiecewiseVanillaModel.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sschlenkrich.github.io/PiecewiseVanillaModel.jl/dev/)
[![Build Status](https://github.com/sschlenkrich/PiecewiseVanillaModel.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sschlenkrich/PiecewiseVanillaModel.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/sschlenkrich/PiecewiseVanillaModel.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sschlenkrich/PiecewiseVanillaModel.jl)


A Piece-wise Model for Vanilla Option Pricing.

The PiecewiseVanillaModel.jl package implements a static model for the pricing and implied volatility calculation of Vanilla options.

The model is based on a piece-wise specification of the terminal distribution of the underlying financial risk factor.

Details of the model are discussed in (TBD).

Model usage is illustrated in the example Julia notebook [JaeckelExample.ipynb](./notebooks/JaeckelExample.ipynb).

An alternative Python implementation with an example is provided in [PiecewiseVanillaModel.ipynb](./python/PiecewiseVanillaModel.ipynb).
