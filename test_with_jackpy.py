import numpy as np

from jackpy.jack import JackPol
from jack import *

def partition_generator(n, max_val=None):
    if max_val is None:
        max_val = n
    if n == 0:
        yield []
        return
    for i in range(min(n, max_val), 0, -1):
        for tail in partition_generator(n - i, i):
            yield [i] + tail

def test_case(deg,part,alpha,
              N=1000,
              low=-2.0,
              high=2.0):
    print(f"Testing {deg=}, {part=}, {alpha=}")

    poly = JackPol(deg,part,alpha=alpha,which='J')
    points = np.random.uniform(low=low,high=high,size=(N,deg))

    sympy_side = np.empty((N,),dtype=np.float64)
    alg_side = np.empty((N,),dtype=np.float64)
    for i in range(N):
        x = poly.eval(list(points[i]))
        y = jack(part,points[i],alpha)
        sympy_side[i] = x
        alg_side[i] = y
    assert np.allclose(sympy_side,alg_side)

def main():
    for deg in range(1,6):
        for alpha in [1.0,1.5,2.0,2.5]:
            for part in partition_generator(deg):
                test_case(deg,part,alpha)

if __name__ == "__main__":
    np.random.seed(3141592)
    main()

