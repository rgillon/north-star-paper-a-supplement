"""Exact certificate for a full pair with binary spectator slots shared by rank."""
from fractions import Fraction
from math import comb, gcd, prod
from sympy import factorint, primerange

from ppi_reallocated_pair import add_damping, donor_bounds
from ppi_reallocated_donors import UPPER, packet_bound
from ppi_reduced_lower import MASKS_1, SIGNED_1, PATH, pinned


def odd_part(n):
    return n // (n & -n)


def blend(bound, target, denominator, dyadic_exponent):
    odd = odd_part(denominator)
    # Scale target by a dyadic multiple of its odd denominator just below one.
    precision = (odd*(1 + int(bound) + 1)).bit_length()
    numerator = (2**precision // odd)*odd
    weight = Fraction(numerator,2**precision)
    assert 0 < weight <= 1 and (1-weight)*bound < 1
    exponent = precision + max(dyadic_exponent, (denominator & -denominator).bit_length()-1)
    return target+1, exponent


def main():
    axes, E, dyadic, upper, lower, U, L = donor_bounds()
    original = pinned(PATH)
    hidden, fibers, masks, scale, determinant = UPPER[0]
    rows = [[i for i in range(len(fibers)) if m >> i & 1] for m in masks]
    d = original["donor"](hidden,fibers,rows,scale,determinant)
    upper, eu = blend(upper,Fraction(packet_bound(original,fibers,rows,d),d),d,0)
    ld,lower_norm = pinned("scripts/ppi_lower_shield_pair.py")["inspect"](SIGNED_1,MASKS_1)
    lower, el = blend(lower,lower_norm,ld,25)
    tags = []
    excluded = prod(axes)
    for p in primerange(3, 100000):
        if gcd(p*(p-1),excluded) == 1:
            tags.append(p)
            excluded *= p
            if len(tags) == 26:
                break
    assert len(tags) == 26
    # Earlier predecessors are smaller than each later tag.
    K = prod(tags)
    aux = prod(odd_part(p-1) for p in tags)
    assert gcd(K*prod(axes),aux) == 1
    r = 13
    val = sorted(((p-1) & -(p-1)).bit_length()-1 for p in tags)
    b, extra = sum(val[:r]), sum(val[-r:])-sum(val[:r])
    M = 2**b
    # Counts by number of nonidentity tag coordinates.
    counts = [1]
    for p in tags:
        counts = [x+y for x,y in zip(counts+[0],[0]+[(p-1)*v for v in counts])]
    n = sum(counts[:r])
    g = (16*K//M).bit_length()
    grid = 2**g
    c = Fraction(((2*n+M)*grid+K-1)//K,grid)
    assert 0 < c < 2
    for t in range(1,5000):
        R = 79+3*t
        s = (16*R*E*K//M).bit_length()
        delta = Fraction(E,2**s)
        error = (4*(1+upper)*(1+lower)+6*R*delta*(1+U)*(1+L))/2**t
        if 8*error < delta:
            break
    else:
        raise ValueError("no verified damping count")
    A = t+max(g+eu+el,s+dyadic)+b+3
    step = Fraction(K,2*M)*delta
    width = Fraction(1)
    while 2*width > step:
        width /= 2
    assert step/4 < width <= step/2
    center = (K*c/2-n)/M
    left = center-R*step-width/2
    right = center+width/2
    assert 0 < left < right < 1
    assert (2**A*left).denominator == 1
    assert (2**A*step).denominator == (2**A*width).denominator == 1
    assert comb(len(tags),r) >= 2*(R+1)*A
    axes.update(tags)
    axes.update(factorint(aux))
    before = len(axes)
    digest = add_damping(axes,t)
    assert len(axes) == before+2*t
    print("tag primes", tags)
    print("auxiliary primes", sorted(factorint(aux)))
    print(f"t={t}, R={R}, g={g}, s={s}, A={A}, extra={extra}, a={A+extra}")
    print(f"slots={comb(len(tags),r)}, demand bound={2*(R+1)*A}")
    print(f"distinct odd primes={len(axes)}, damping SHA-256={digest}")


if __name__ == "__main__":
    main()
