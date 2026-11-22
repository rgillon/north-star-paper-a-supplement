"""Exact partial upper interpolation by adding two balanced singleton blocks."""
from fractions import Fraction
from math import isqrt, prod
from sympy import factorint, primerange

from ppi_mixed_slots import verify_templates
from ppi_rank_box_pair import odd_part
from ppi_weighted_fibers import HIDDEN, donor_bounds, tag_parameters


FIBERS = [(23, 5, 1, 18446462601427681280), (-31, 2, 2, 11565243845243305984),
 (37, 2, 2, 12144106481226088448), (-41, 2, 1, 16140901066643341312),
 (-43, 2, 1, 14467813805217808384), (-47, 2, 1, 16204127933579689984),
 (67, 1, 1, 11565420317388079104), (71, 1, 1, 12249931726754676736),
 (79, 1, 1, 12297264190501980160), (101, 1, 1, 18374686482355978240),
 (-107, 1, 1, 18412967080698576896), (-109, 1, 1, 18437736877147553792),
 (-113, 1, 1, 18442381215611748352)]
MASKS = [1037, 600, 4614, 816, 3585, 5444, 2251, 630, 4435, 1641, 1713, 1930]
PARTIAL = ((43781, 43783, 43787, 43789), (457, 461, 463, 467))


def apply_partial_blocks(base, blocks):
    axes,E,dy,upper,lower,U,L,eu,el,fixed,pad = base
    for primes in blocks:
        assert len(primes) in (2,4) and list(primes) == sorted(set(primes))
        if len(primes) == 4:
            a,b,c,d = primes
            assert a+d == b+c
        assert axes.isdisjoint(primes)
        assert all(p > 2 and all(p % i for i in range(2,isqrt(p)+1)) for p in primes)
        Z = prod(q-1 for q in primes)
        # The two scales minimize the endpoint spread and balance the root
        # cost against the old norm times the largest singleton deficit.
        low,high = primes[0]-1,primes[-1]-1
        scale = min(Fraction(2*low*high,low+high), U*high/(U+high))
        bits = Z.bit_length()+24
        alpha = Fraction((1 << bits)*scale//Z, 1 << bits)
        root = alpha*Z
        error = Fraction(0)
        for mask in range(1,1 << len(primes)):
            packet = prod(1 if mask >> i & 1 else 1-q for i,q in enumerate(primes))
            error = max(error,abs(int(mask.bit_count()==1)+alpha*packet))
        U = max(Fraction(1),root,U*error)
        axes.update(primes)
        dy += bits
        fixed += len(primes)
    # The balanced-block midpoint anchor spends one bit per block. This
    # also bounds the support-indicator extension for unbalanced pairs.
    return axes,E,dy,upper,lower,U,L,eu+len(blocks),el,fixed,pad


def partial_bounds():
    base = donor_bounds(HIDDEN,FIBERS,MASKS,8832,113,73572 << 67,(78,1340,64,1354))
    assert base[5] < 428778222
    result = apply_partial_blocks(base,PARTIAL)
    assert result[5] < 475
    return result


def main():
    args = partial_bounds()
    raw = [p for p in primerange(3,1000) if p not in args[0]][:9]
    safe = [p for p in primerange(3,30000) if p not in args[0]|set(raw)
            and set(factorint(odd_part(p-1))) <= set(raw)]
    data = tag_parameters(args,raw,safe[4:18],safe[:4],7,3)
    total = verify_templates(args,*data,650,33,Fraction(8,3),head_count=252)
    assert total == 142


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
