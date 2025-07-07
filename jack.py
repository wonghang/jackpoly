# This code is translated from jack.m written by Plamen Koev 2003 with the help of Grok 3
#
# function f=jack(l,x,alpha). Computes the Schur function corresponding
# to partition l with variables specified in x.
# Plamen Koev, 2003
#
import numpy as np
from typing import Sequence

__all__ = [
    'jack',
]

def jack(l: Sequence[int],
         xx: Sequence[float],
         alpha1: float) -> float:
    """
    Evaluate the J-normalized Jack polynomial of a given partition and alpha parameter.

    - `l` is a partition represented by a list of positive integers.
    - `alpha1` is the alpha parameter.
    - `xx` is a list or NumPy one-dimensional array of evaluation points.

    Parameters
    ----------
    l : array_like of int
        A partition represented as a non-increasing sequence of positive integers
        (e.g., [3, 2, 2] for the partition 3 + 2 + 2).
        Zeros or negative values in the list are ignored.

    xx : array_like of float
        A one-dimensional array or iterable of real numbers representing the variables
        at which the Jack polynomial is evaluated. Must be convertible to a 1D NumPy array.

    alpha1 : float
        The alpha parameter

    Returns
    -------
    float
        The evaluated value of the J-normalized Jack polynomial at the point

    Raises
    ------
    ValueError
        If `xx` is not one-dimensional or not a floating-point array.
    """
    x = np.asarray(xx,dtype=np.float64)

    if x.ndim != 1:
        raise ValueError("xx should be of one-dimensional")

    alpha = np.float64(alpha1)
    n = len(x)

    l = np.array(l,dtype=int)
    l = l[l > 0]
    Lmax = l + 1
    Lp = len(l)

    ja = np.full((np.prod(Lmax) + 1, n), np.inf, dtype=np.float64)

    lma = np.zeros(Lp, dtype=int)
    if Lp > 0:
        lma[Lp-1] = 1
        for i in range(Lp-2, -1, -1):
            lma[i] = lma[i+1] * Lmax[i+1]

        result = jack1(n, 0, nmu(l, Lmax, Lp), nmu(l, Lmax, Lp), x, Lp, Lmax, n, ja, lma, alpha)
    else:
        result = np.float64(1.0)

    return result

def part(pn, i, lma):
    """
    given an integer pn that represents a partition, this function returns
    part i of the partition in O(1) time
    """
    if i > len(lma):
        return 0
    if i == 1:
        return pn // lma[i-1]
    else:
        return (pn % lma[i-2]) // lma[i-1]

def nmu(l, Lmax, Lp):
    """
    nmu computes the unique integer that represents the partition l
    """
    f = 0
    for i in range(Lp):
        f = Lmax[i] * f
        if i < len(l):
            f += l[i]
    return f

def jack1(j, k, lambda_, l, x, Lp, Lmax, n, ja, lma, alpha):
    """
    Given vector x and partition, represented by an integer l, jack1(j,k,lambda,l)
    computes the Jack polynomial J_lambda^alpha (x_1,...,x_j), where the partitions
    l, such that lambda-l is horizontal strip, are generated to have
    lambda(i)=l(i) for i=1,2,...,k. The parameter k is designed to keep
    track of the recursion. from outside the function should be called
    with ja(j,0,lambda,lambda). This function can only be called after an initialization
    by the function jack
    """
    s = 1.0

    if 1 <= j and l > 0:
        t = ja[l, j-1]
        if k == 0 and t != np.inf:
            s = t
        elif part(l, j+1, lma) > 0:
            s = 0.0
        elif j == 1:
            p = part(l, 1, lma)
            s = x[0]**p * np.prod(1 + alpha * np.arange(p))
        else:
            i = k if k > 0 else 1
            s = jack1(j-1, 0, l, l, x, Lp, Lmax, n, ja, lma, alpha) * AB(lambda_, l, Lp, lma, alpha) * x[j-1]**lm(lambda_, l, Lp, lma)
            while part(l, i, lma) > 0:
                if part(l, i, lma) > part(l, i+1, lma):
                    if part(l, i, lma) > 1:
                        s += jack1(j, i, lambda_, l - lma[i-1], x, Lp, Lmax, n, ja, lma, alpha)
                    else:
                        s += (jack1(j-1, 0, l - lma[i-1], l - lma[i-1], x, Lp, Lmax, n, ja, lma, alpha) *
                              AB(lambda_, l - lma[i-1], Lp, lma, alpha) * x[j-1]**lm(lambda_, l - lma[i-1], Lp, lma))
                i += 1
        if k == 0:
            ja[l, j-1] = s
    return s

def AB(l, m, Lp, lma, alpha):
    """
    Computes the coefficient A_lambda/B_mu for Jack polynomials.
    """
    f = alpha.dtype.type(1.0)
    for i in range(1, Lp+1):
        for j in range(1, part(m, i, lma)+1):
            if l_t(l, j, lma) == l_t(m, j, lma):
                f /= (l_t(m, j, lma) - i + alpha * (part(m, i, lma) - j + 1))
            else:
                f /= (l_t(m, j, lma) - i + 1 + alpha * (part(m, i, lma) - j))
    for i in range(1, Lp+1):
        for j in range(1, part(l, i, lma)+1):
            if l_t(l, j, lma) == l_t(m, j, lma):
                f *= (l_t(l, j, lma) - i + alpha * (part(l, i, lma) - j + 1))
            else:
                f *= (l_t(l, j, lma) - i + 1 + alpha * (part(l, i, lma) - j))
    return f

def lm(l, m, Lp, lma):
    """
    computes |lambda/mu|
    """
    f = 0
    for i in range(1, Lp+1):
        f += part(l, i, lma) - part(m, i, lma)
    return f

def l_t(l, q, lma):
    """
    compute the q-th part of the transpose of a partition lt(l,q)=#(l(i)|l(i)>=q)
    """
    i = 1
    f = 0
    while part(l, i, lma) >= q:
        f += 1
        i += 1
    return f

