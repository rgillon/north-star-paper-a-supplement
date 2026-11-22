"""Exact weighted upper donor and compressed Level 1 pair."""
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, lcm, prod
from sympy import Matrix, factorint, primerange

from ppi_mixed_slots import verify_templates, weight_counts
from ppi_rank_box_pair import odd_part
from ppi_modular_pair import assemble_bounds

HIDDEN = (3, 5, 7, 11, 17, 53)
FIBERS = [(23, 7, 18341059568487923712), (-31, 5, 18374932773497470976),
 (-37, 3, 12294968272247324672), (-41, 4, 18445195963619344384),
 (43, 3, 16177106335278563328), (-47, 3, 17357013704067022848),
 (-67, 2, 16204127934663917568), (71, 2, 17870468043520376832),
 (73, 2, 18059434508045778944), (79, 2, 18442240478322622464),
 (-83, 2, 18446708892030238720), (101, 1, 12249940522176610304),
 (107, 1, 12144106481369254016), (-109, 1, 12294968272247324672),
 (113, 1, 12297829382470805632), (-127, 1, 16177106335278563328),
 (-131, 1, 16195191100384739328), (137, 1, 17348006504275378176),
 (139, 1, 17357013704067022848), (149, 1, 18338657685504786432)]
MASKS = [523, 284, 1313, 131094, 131145, 368, 36897, 263233, 81940, 278564, 9488, 1728, 294960,
 525856, 10432, 557192, 197128, 165892, 23554]
SCALE = 6192
MINOR = -298
LOWER = [[59, 61, 163, 173, 223, 233, 281, 293, -151, -199, -229, -277, -311, -317],
 [3095, 4049, 5687, 6455, 6743, 8731, 8743, 8848, 9044, 12333, 13427, 13511, 15486]]


def packet_extrema(primes, masks, coefficients, lower=False):
    """Integer numerator at each order, by the tensor product of local evaluations."""
    values = [0]*(1 << len(primes))
    values[0] = -sum(coefficients)
    for mask, coefficient in zip(masks,coefficients):
        values[mask] += coefficient
    for i,p in enumerate(primes):
        stride = 1 << i
        for start in range(0,len(values),2*stride):
            for j in range(start,start+stride):
                root, nontrivial = values[j], values[j+stride]
                values[j] = root+nontrivial
                if lower:
                    quotient, remainder = divmod(nontrivial,p-1)
                    assert remainder == 0
                    values[j+stride] = root-quotient
                else:
                    values[j+stride] = root+(1-p)*nontrivial
    return min(values),max(values)


def solve_packets(n, masks, rhs):
    matrix = Matrix([[m >> i & 1 for m in masks] for i in range(n)])
    coefficients = matrix[:-1,:].inv()*Matrix(rhs[:-1])
    assert matrix*coefficients == Matrix(rhs)
    return [Fraction(a) for a in coefficients]


def donor_bounds(hidden, data, conductors, scale, minor, upper_multiplier, lower_offsets,
                 hidden_powers=None):
    primes = [abs(p) for p,_,_,_ in data]
    n = len(primes)
    assert all(p > 2 and all(p%i for i in range(2,isqrt(p)+1)) for p in list(hidden)+primes)
    assert all(c > 0 for _,c,_,_ in data)
    assert sum(c if p > 0 else -c for p,c,_,_ in data) == 0
    powers = [1]*len(hidden) if hidden_powers is None else hidden_powers
    assert len(powers) == len(hidden) and all(e >= 1 for e in powers)
    strides = [prod(e+1 for e in powers[:i]) for i in range(len(hidden))]
    states = range(prod(e+1 for e in powers))
    levels = [[s//d % (e+1) for d,e in zip(strides,powers)] for s in states]
    masses = [prod((p-1)*p**(e-1) for p,e in zip(hidden,row) if e) for row in levels]
    for p,c,exponent,mask in data:
        assert exponent in (1,2)
        assert 0 < mask < 1 << len(masses)
        shells = [s for s in states if mask >> s & 1]
        ratio = abs(p)-1 if exponent == 1 else abs(p)
        assert sum(masses[s] for s in shells) == scale*c*ratio
        assert all(mask >> (s+d) & 1 for s in shells for i,d in enumerate(strides)
                   if levels[s][i] < powers[i])
    assert len(conductors) == n-1 and all(0 < m < 1 << n for m in conductors)
    rows = [[i for i in range(n) if mask >> i & 1] for mask in conductors]
    assert all(sum(data[i][0]*data[i][1] for i in row) == 0 for row in rows)
    matrix = [[int(i in r) for r in rows] for i in range(n-1)]
    assert Matrix(matrix).det() == minor
    du = abs(minor)*prod(primes)
    fibers = [(p,mask) for p,_,_,mask in data]
    co = solve_packets(n,conductors,[-du//p for p in primes])
    assert all(v.denominator == 1 for v in co)
    ub = Fraction(max(map(abs,packet_extrema(primes,conductors,list(map(int,co))))),du)
    signed,masks = LOWER
    ps = list(map(abs,signed))
    assert all(p > 2 and all(p%i for i in range(2,isqrt(p)+1)) for p in ps)
    assert sum((p-1)*(s//p) for p,s in zip(ps,signed)) == 0
    assert len(masks) == len(ps)-1 and all(0 < m < 1 << len(ps) for m in masks)
    assert all(sum(s for i,s in enumerate(signed) if m >> i & 1) == 0 for m in masks)
    coeff = solve_packets(len(ps),masks,[Fraction(1-p,p) for p in ps])
    dl = lcm(sum(coeff).denominator,
             *(v.denominator*prod(p-1 for i,p in enumerate(ps) if m >> i & 1)
               for v,m in zip(coeff,masks)))
    scaled = [dl*a for a in coeff]
    assert all(v.denominator == 1 for v in scaled)
    lb = Fraction(max(map(abs,packet_extrema(ps,masks,list(map(int,scaled)),True))),dl)
    return assemble_bounds(hidden,fibers,signed,du,ub,dl,lb,upper_multiplier,lower_offsets)


def bounds():
    data = [(p,c,1,mask) for p,c,mask in FIBERS]
    result = donor_bounds(HIDDEN,data,MASKS,SCALE,MINOR,72856,(64,2358,0,2422))
    assert result[5] < 1746300 and result[6] < Fraction(2447,1000)
    return result


def tag_parameters(args, raw, fine, common, r, d):
    tags = raw + common + fine
    assert len(tags) == len(set(tags)) and args[0].isdisjoint(tags)
    assert all(p > 2 and all(p%i for i in range(2,isqrt(p)+1)) for p in tags)
    assert len(common)+r < len(tags)-d
    aux = prod(odd_part(p-1) for p in common+fine)
    assert gcd(aux,prod(args[0])*prod(common+fine)) == 1
    full = prod(p-1 for p in tags)
    coarse = sorted((full//prod(p-1 for p in omitted) for omitted in combinations(tags,d)
                     if any(p in common for p in omitted)),reverse=True)
    fixed_mass = prod(p-1 for p in common)
    fine_weights = [fixed_mass*prod(p-1 for p in slot) for slot in combinations(fine,r)]
    k = len(tags)-d
    n = (sum(weight_counts(tags)[:k])-fixed_mass*sum(weight_counts(raw+fine)[:k-len(common)])
         +fixed_mass*sum(weight_counts(fine)[:r]))
    return tags,aux,prod(tags),n,coarse,fine_weights


def main():
    args = bounds()
    # Even without a spectator cost, the present base-error envelope cannot
    # fit the 100-prime budget with these fixed blocks and two primes per pad.
    U,L = args[5:7]
    base = (1+U)*(3+L)/4
    assert all(base*(args[-2]+2*t) > 2**(t+1)
               for t in range((100-len(args[0]))//2+1))
    raw = [p for p in primerange(3,1000) if p not in args[0]][:7]
    safe = [p for p in primerange(3,30000) if p not in args[0]|set(raw)
            and set(factorint(odd_part(p-1))) <= set(raw)]
    common, fine = safe[:4], safe[4:18]
    data = tag_parameters(args,raw,fine,common,7,3)
    total = verify_templates(args,*data,700,48,Fraction(16,5),head_count=320)
    assert total == 169


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
