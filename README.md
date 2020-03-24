# DiskArrayTools.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
[![Build Status](https://travis-ci.com/meggart/DiskArrayTools.jl.svg?branch=master)](https://travis-ci.com/meggart/DiskArrayTools.jl)
[![codecov.io](http://codecov.io/github/meggart/DiskArrayTools.jl/coverage.svg?branch=master)](http://codecov.io/github/meggart/DiskArrayTools.jl?branch=master)

This package collects additional functionalities to work with DiskArrays.
Currently it provides:

- `diskstack` concatenating arrays of the same type to a new DiskArray along a new dimension
- `InterpolatedDiskArray` interpolates data to a new grid using Interpolations.jl
