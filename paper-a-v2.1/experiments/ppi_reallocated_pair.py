"""Exact parameters and finite block checks for the compressed full pair."""
from fractions import Fraction
from hashlib import sha256
from itertools import islice
from math import comb, gcd, prod
from sympy import isprime, primerange

from ppi_reallocated_donors import UPPER, packet_bound
from ppi_reduced_lower import SIGNED_1, SIGNED_2, MASKS_1, MASKS_2, pinned, PATH


def donor_bounds():
    upper = pinned(PATH)
    lower = pinned("scripts/ppi_lower_shield_pair.py")
    axes = set(map(abs,SIGNED_1 + SIGNED_2))
    du, individual = [], []
    for hidden, fibers, masks, scale, det in UPPER:
        rows = [[i for i in range(len(fibers)) if m >> i & 1] for m in masks]
        d = upper["donor"](hidden,fibers,rows,scale,det)
        du.append(d)
        individual.append(Fraction(packet_bound(upper,fibers,rows,d),d))
        axes.update(hidden)
        axes.update(abs(p) for p,_ in fibers)
    assert len(axes) == 79 and gcd(*du) == 1
    k = du[0]*pow(du[0],-1,du[1])
    if 2*k > prod(du):
        k -= prod(du)
    upper_bound = abs(k)*individual[0] + abs(1-k)*individual[1]
    (d1,l1),(d2,l2) = [lower["inspect"](c,m) for c,m in
                        [(SIGNED_1,MASKS_1),(SIGNED_2,MASKS_2)]]
    u,v = d1//(d1&-d1),d2//(d2&-d2)
    assert gcd(u,v) == 1
    assert max(d1 & -d1,d2 & -d2) <= 2**25
    k = u*pow(u,-1,v)
    lower_bound = abs(k)*l1 + abs(1-k)*l2
    E = prod(du)*d1*d2
    E //= E & -E
    dyadic = max((d & -d).bit_length()-1 for d in du) + max(
        (d & -d).bit_length()-1 for d in (d1,d2))
    return axes, E, dyadic, upper_bound, lower_bound, max(individual), max(l1,l2)


def add_damping(axes, t):
    digest = sha256()
    count = 0
    for p in primerange(2**21,2**24):
        q = p*(p-1)+1
        if p in axes or q in axes or not isprime(q):
            continue
        assert q < 2**64
        axes.update((p,q))
        digest.update(f"{p},{q}\n".encode("ascii"))
        count += 1
        if count == t:
            break
    assert count == t
    return digest.hexdigest()


def parameters():
    axes, E, dyadic, upper_bound, lower_bound, U, L = donor_bounds()
    h, m = 160, 2**18
    pairs = ((p,m+2-p) for p in primerange(3,(m+2)//2)
             if p not in axes and m+2-p not in axes and isprime(m+2-p))
    tags = [p for pair in islice(pairs,2*h) for p in pair]
    assert len(set(tags)) == 4*h
    K, M = prod(tags), m**h
    background = (1+2*m)**h-(2*m)**h
    g = (16*K//M).bit_length()
    grid = 2**g
    c = Fraction(((2*background+M)*grid+K-1)//K,grid)
    assert 0 <= c-Fraction(2*background+M,K) < Fraction(M,16*K)
    for t in range(5000,10001):
        R = 79+3*t
        s = (16*R*E*K//M).bit_length()
        delta = Fraction(E,2**s)
        err_times_power = (4*(1+upper_bound)*(1+lower_bound)
                           +6*R*delta*(1+U)*(1+L))
        if 3*err_times_power >= delta*2**t:
            continue
        a = t+max(g+25,s+dyadic)+18*h+1
        if 2**79*comb(t+5,5)*a <= 2**h:
            break
    else:
        raise ValueError("parameter search exhausted")
    assert Fraction(M,32*R*K) <= delta <= Fraction(M,16*R*K)
    assert 3*err_times_power < delta*2**t
    assert 2*background < K*(c-delta*R-delta/3)
    assert K*(c+delta/3) < 2*(background+M)
    axes.update(tags)
    digest = add_damping(axes,t)
    assert len(axes) == 79+4*h+2*t
    print(f"h={h}, t={t}, a={a}, g={g}, s={s}, R={R}, distinct odd primes={len(axes)}")
    print("Damping list SHA-256",digest)


if __name__ == "__main__":
    parameters()
