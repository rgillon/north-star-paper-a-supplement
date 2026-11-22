"""Reproduce both prime tuples, exact helper witnesses and the full pair."""
from fractions import Fraction
from math import gcd, isqrt, prod
from sympy import factorint, primerange

from ppi_block_rank_pair import UPPER_DONORS
from ppi_coarse_pair import select_sections
from ppi_rank_box_pair import blend, odd_part
from ppi_reallocated_donors import packet_bound
from ppi_reduced_lower import MASKS_1, SIGNED_1, PATH, pinned


def pocklington(n, known_primes):
    F = 1
    for q in sorted(set(known_primes)):
        assert (n-1)%q == 0
        m = n-1
        while m%q == 0:
            F *= q
            m //= q
        assert any(pow(a,n-1,n)==1 and gcd(pow(a,(n-1)//q,n)-1,n)==1
                   for a in range(2,256))
    assert F > isqrt(n)


def modular_bounds():
    h,f,masks,scale,det = UPPER_DONORS[0]
    routine = pinned(PATH)
    rows = [[i for i in range(len(f)) if m >> i & 1] for m in masks]
    du = routine["donor"](h,f,rows,scale,det)
    ub = Fraction(packet_bound(routine,f,rows,du),du)
    dl,lb = pinned("scripts/ppi_lower_shield_pair.py")["inspect"](SIGNED_1,MASKS_1)
    return assemble_bounds(h,f,SIGNED_1,du,ub,dl,lb,76751,(391,547,91,847))


def assemble_bounds(h,f,signed_lower,du,ub,dl,lb,upper_multiplier,lower_offsets):
    ol = odd_part(dl)
    k = upper_multiplier*prod(primerange(2,32))*du
    upper = [1+k*a for a in (2,12,3,4)]
    small = list(map(int,factorint(12*k)))
    assert all(2 <= q <= 100000 and all(q%i for i in range(2,isqrt(q)+1)) for q in small)
    for q in upper:
        pocklington(q,[p for p in small if (q-1)%p==0])
    P = upper[0]
    lower = [1+P*(pow(P,-1,ol)+j*ol) for j in lower_offsets]
    for q in lower:
        small_factors = [p for p in primerange(2,1000) if (q-1)%p == 0]
        pocklington(q,[P]+small_factors)
    axes = list(h)+[abs(p) for p,_ in f]+list(map(abs,signed_lower))+upper+lower
    assert len(set(axes)) == len(axes) and sorted(axes)[-4:] == sorted(lower)
    signs = (1,1,-1,-1)
    assert sum(Fraction(s,q-1) for s,q in zip(signs,upper)) == 0
    assert sum(s*(q-1) for s,q in zip(signs,lower)) == 0
    bounds = []
    for qs,modulus,is_upper in [(upper,du,True),(lower,ol,False)]:
        norm_v = norm_c = 0
        for active in range(16):
            packet = prod(1 if active >> i & 1 else 1-q for i,q in enumerate(qs))
            if is_upper:
                v = 1-packet
                forced = active.bit_count() >= 3
                target = int(active.bit_count() == 3)
            else:
                v = pow(-2,-1,ol)*(packet-prod(1-q for q in qs))
                forced = active.bit_count() <= 1
                target = int(active.bit_count() == 1)
            norm_v = max(norm_v,abs(v))
            if forced:
                assert (target-v)%modulus == 0
                norm_c = max(norm_c,abs((target-v)//modulus))
        bounds.append((norm_v,norm_c))
    v,c = bounds[0]
    U,eu = blend(du*ub*c+v,ub,du,0)
    v,c = bounds[1]
    L,el = blend(ol*lb*c+v,lb,dl,(dl&-dl).bit_length()-1)
    fixed_rank = len(f)+2+len(signed_lower)+8
    E = odd_part(du*dl)
    dyadic = ((du&-du)*(dl&-dl)).bit_length()-1
    return set(axes),E,dyadic,U,L,ub,lb,eu,el,fixed_rank,2


def main():
    args = modular_bounds()
    print("fixed axes",len(args[0]),"fixed rank bound",args[-2])
    select_sections(*args)


if __name__ == "__main__":
    main()
