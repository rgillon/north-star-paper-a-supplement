"""Exact five-prime example for the private marginal and exposed shared shells."""
from fractions import Fraction
from itertools import product
from math import prod


PRIVATE = (11, 13)
SHARED = (3, 5, 7)
# Shared order exponents are one plus these triples. Each list belongs to
# the indicated private activity mask, low bit first.
DIGITS = {
    0: [(2,3,2), (5,1,1), (1,0,6), (0,4,3), (1,5,1),
        (3,1,3), (1,3,3), (0,6,1), (0,2,5), (4,0,3)],
    1: [(5,2,0), (3,2,2), (3,4,0), (4,1,2), (7,0,0), (3,0,4)],
    2: [(2,2,3), (4,3,0), (2,5,0), (6,1,0)],
    3: [(5,0,2)],
}


def shell_size(prime, exponent):
    return (prime-1)*prime**(exponent-1) if exponent else 1


def shell_sum(prime, bound, exponent, conductor):
    """Sum over the spatial shell, by subtracting its two subgroup sums."""
    def subgroup(v):
        return prime**v if conductor <= bound-v else 0
    return subgroup(exponent) - (subgroup(exponent-1) if exponent else 0)


def main():
    primes = PRIVATE+SHARED
    labels = {}
    for mask, triples in DIGITS.items():
        for triple in triples:
            label = tuple(v+1 for v in triple)
            assert label not in labels and sum(label) == 10
            labels[label] = tuple((mask >> i) & 1 for i in range(2))
    bounds = (1,1)+tuple(map(max,zip(*labels)))
    support = {private+label: (-1)**sum(private) for label,private in labels.items()}
    orders = list(product(*(range(a+1) for a in bounds)))
    def below(v, w):
        return all(i <= j for i,j in zip(v,w))
    strict_lower = {v for v in orders if v not in support and any(below(v,w) for w in support)}
    # Check both order ideals directly on the full divisor lattice.
    for sign in (1,-1):
        ideal = strict_lower | {v for v,s in support.items() if s == sign}
        for v in ideal:
            for j,a in enumerate(v):
                if a:
                    assert v[:j]+(a-1,)+v[j+1:] in ideal

    def transform(conductor):
        return sum(s*prod(shell_sum(p,a,v,k) for p,a,v,k in zip(primes,bounds,w,conductor))
                   for w,s in support.items())
    mixed = {k:transform(k) for k in orders}
    private_rows = [mixed[k+(0,0,0)] for k in product(range(2),repeat=2)]
    assert private_rows == [0,0,0,11907*48*11*13]
    masses = [sum(prod(shell_size(p,v) for p,v in zip(primes,w))
                  for w in support if w[:2] == a) for a in product(range(2),repeat=2)]
    assert masses == [11907*48*120]*4

    for e,x in labels.items():
        h = tuple(a-b+1 for a,b in zip(bounds[2:],e))
        assert all(mixed[k+h] != 0 for k in product(range(2),repeat=2))
        weight = {}
        for y in product(*(range(b+1) for b in e)):
            weight[y] = prod(Fraction(shell_sum(p,a,v,k),2*p**(b-1))
                             for p,a,v,k,b in zip(SHARED,bounds[2:],y,h,e))
        assert sum(weight.values()) == 0 and sum(map(abs,weight.values())) == 1
        # The character packet has constant value on each spatial shell.
        support_moment = sum(w for y,w in weight.items() if x+y in support)
        upper_moment = sum(w for y,w in weight.items() if any(below(v,x+y) for v in support))
        assert support_moment == upper_moment == Fraction(-1,8)
        forced = {y:1 if x+y in support else 2 if x+y in strict_lower else None for y in weight}
        assert None not in forced.values()
        assert sum(w*forced[y] for y,w in weight.items()) == Fraction(1,8)

    print(f"{len(support)} support shells, {len(mixed)} conductor tuples, private shell mass {masses[0]}.")
    print(f"Private-only Fourier rows {private_rows}.")
    print("Every exposed shared shell gives support and upper moments -1/8, forced profile moment 1/8.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
