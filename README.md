# lagrange-interpolation

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

Simple library for Lagrange Interpolation in 1D, 2D, and 3D on a
Cartesian product grid. Written in C for maximum compatibility. The
plot below shows the error for piecewise 2nd-order using this library:

![2nd-order interpolation](hthttps://raw.githubusercontent.com/Yurlungur/lagrange-interpolation/master/images/test_fo.png)

## Installation

This code is simple enough that you can just copy-paste `lagrange.c`
and `lagrange.h` into your project. At the moment, that's the only
installation method.

## Use

The header file should provide the interface for the functions you are
interested in. Most interesting, probably, are the
`lagrange_interp_NDfo` functions, where you specify an order of
interpolating polynomial in each dimension and the library evaluates a
piecewise interpolating polynomial of that order over the grid. The
`lagrange_interp_ND` series of functions are global interpolators and
should be used only if your grid points are stable for high-order
interpolation.

## Testing

You can test the code by cloning the directory, entering it, and
typing `make test`. You can also generate a nice plot with `make
plot`. The plotting functionality requires the scientific python
stack.

```bash
git cone git@github.com:Yurlungur/lagrange-interpolation.git
cd lagrange-interpolation
make test
```
