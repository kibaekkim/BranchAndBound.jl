# BranchAndBound.jl
[![Build Status](https://travis-ci.com/kibaekkim/BranchAndBound.jl.svg?token=3N6HLyM8rqygf5Rmoqzp&branch=master)](https://travis-ci.com/kibaekkim/BranchAndBound.jl)
[![codecov](https://codecov.io/gh/kibaekkim/BranchAndBound.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kibaekkim/BranchAndBound.jl)

This Julia package provides an abstract framework for branch-and-bound methods. Applications of this package include:

- a branch-and-bound method for mixed-integer linear programming (see `./examples/milp.jl`)
- a branch-and-bound method for mixed-integer nonlinear programming (see `./examples/minlp.jl`)


## Installation

This package can be installed by

```julia
] add https://github.com/kibaekkim/BranchAndBound.jl.git
```

## Acknowledgement

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.
