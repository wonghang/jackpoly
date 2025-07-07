### Introduction

This is a Python implementation for efficiently computing the (J-normalized) [Jack polynomial](https://mathworld.wolfram.com/JackPolynomial.html), based on the algorithm described in the paper:

*Koev, Plamen, and Alan Edelman. "The efficient evaluation of the hypergeometric function of a matrix argument." Mathematics of Computation 75.254 (2006): 833–846.*

The original code comes from *jack.m* by Prof. Plamen Koev, on his [website](https://sites.google.com/sjsu.edu/plamenkoev/home/software).

I translated it into Python + NumPy with the help of [Grok](https://x.ai/grok).

### Prerequisites

[numpy](https://numpy.org/)>=1.23.3

[jackpy](https://github.com/stla/jackpy) (if you want to run `test_with_jackpy.py`)

### License

GNU General Public License v3.0
