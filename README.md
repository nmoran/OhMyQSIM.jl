# OhMyQSIM

![](https://travis-ci.com/nmoran/OhMyQSIM.jl.svg?token=CWy6KXR9ECn794Hyhkpx&branch=master)

Yet another quantum simulator/learning experience.

## Installation and testing
To install the package to your current environment use:

```
using Pkg
Pkg.add(PackageSpec(url="https://github.com/nmoran/OhMyQSIM.jl"))
```

The unittests can be run using
```
Pkg.test("OhMyQSIM")
```

## Basic functionality

A simple example is as follows

```
using OhMyQSIM

qreg = FullStateQuantumRegister{ComplexF64}(3, "000")
apply_1qubit!(qreg, Gates.x, 1)
to_str(qreg)
```

which gives

```
"(1.0 + 0.0im|100>)"
```
