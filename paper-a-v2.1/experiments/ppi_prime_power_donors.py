"""Exact lower interpolation with eight axes and the coatom reciprocal bound."""
from fractions import Fraction
from itertools import product
from math import isqrt, lcm, prod
from sympy import Matrix


PRIMES = (3,5,7,13,17,19,23,53)
POWERS = (2,2,2,1,1,1,1,1)
SIGNS = (1,1,-1,1,1,1,1,-1)
MASKS = (60,85,163,216,836,1328,1352,1545,1042,1120)


def verify_lower():
    n = len(PRIMES)
    mass = [p**(a-1)*(p-1) for p,a in zip(PRIMES,POWERS)]
    assert sum(s*m for s,m in zip(SIGNS,mass)) == 0
    conductors = []
    for mask in MASKS:
        assert 0 < mask < 1 << (n+3)
        levels = [mask >> i & 1 for i in range(n)]
        for i in range(3):
            if mask >> (n+i) & 1:
                assert levels[i] == 0
                levels[i] = 2
        drop = [0 if k == 0 else p**a if k == 1 else m
                for p,a,m,k in zip(PRIMES,POWERS,mass,levels)]
        assert sum(s*v for s,v in zip(SIGNS,drop)) == 0
        conductors.append(levels)
    incidence = Matrix([[mask >> i & 1 for mask in MASKS] for i in range(n+3)])
    rhs = Matrix([Fraction(1-p,p) for p in PRIMES]+[0]*3)
    rows = list(incidence.T.rref()[1])
    assert len(rows) == len(MASKS)
    co = incidence[rows,:].inv()*rhs[rows,:]
    assert incidence*co == rhs
    co = list(map(Fraction,co))
    norm,den = 0,1
    for order in product(*(range(a+1) for a in POWERS)):
        value = -sum(co)
        for alpha,levels in zip(co,conductors):
            packet = prod(1 if v+k <= a else Fraction(-1,p-1) if v+k == a+1 else 0
                          for v,k,p,a in zip(order,levels,PRIMES,POWERS))
            value += alpha*packet
        if sum(bool(v) for v in order) <= 1:
            assert value == int(any(v == a for v,a in zip(order,POWERS)))
        norm = max(norm,abs(value))
        den = lcm(den,value.denominator)
    assert norm < Fraction(241,100)
    return norm,den


def main():
    norm,den = verify_lower()
    primes = [p for p in range(3,1224,2) if all(p % d for d in range(2,isqrt(p)+1))]
    assert len(primes) == 199
    assert sum(Fraction(1,p-1) for p in primes[:-1]) < 2
    assert sum(Fraction(1,p-1) for p in primes) > 2
    print(f"Eight-prime witness, 864 order shells, norm={norm}, denominator={den}.")
    print("The predecessor reciprocal sum first exceeds two at the 199th odd prime.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
