"""Prime-power spectator slots with binary decoding of odd digit weights."""
from itertools import combinations_with_replacement
from fractions import Fraction
from math import gcd, prod
from sympy import factorint, primerange

from ppi_mixed_slots import verify_templates
from ppi_partial_upper_pair import partial_bounds
from ppi_rank_box_pair import odd_part


def shell(p, exponent):
    return (p-1)*p**(exponent-1) if exponent else 1


def rank_masses(primes, bounds, grades):
    masses = [1]
    for p, bound, grade in zip(primes,bounds,grades):
        next_masses = [0]*(len(masses)+bound*grade)
        for degree, value in enumerate(masses):
            for power in range(bound+1):
                next_masses[degree+grade*power] += value*shell(p,power)
        masses = next_masses
    return masses


def same_grade(bounds, grades, target):
    if not bounds:
        if target == 0:
            yield ()
        return
    for e in range(min(bounds[0],target//grades[0])+1):
        for rest in same_grade(bounds[1:],grades[1:],target-e*grades[0]):
            yield (e,)+rest


def spectator_data(axes, raw, guards, fine, raw_bounds, fine_grade, rank, deficit):
    tags = raw+guards+fine
    bounds = raw_bounds+[1]*(len(guards)+len(fine))
    K = prod(p**e for p,e in zip(tags,bounds))
    assert len(tags) == len(set(tags)) and axes.isdisjoint(tags)
    # B contains the full raw-prime part of every tag predecessor product.
    predecessors = prod(p-1 for p in tags)
    factors = factorint(predecessors)
    aux = prod(p**factors.get(p,0) for p in raw)
    assert gcd(aux,prod(guards+fine)) == 1
    guard_mass = prod(p-1 for p in guards)
    cut = sum(bounds)-deficit
    sub_bounds = raw_bounds+[1]*len(fine)
    grades = [1]*len(raw)+[fine_grade]*len(fine)
    assert len(guards)+rank < cut
    background = (sum(rank_masses(tags,bounds,[1]*len(tags))[:cut])
                  -guard_mass*sum(rank_masses(raw+fine,sub_bounds,[1]*len(sub_bounds))[:cut-len(guards)])
                  +guard_mass*sum(rank_masses(raw+fine,sub_bounds,grades)[:rank]))
    ordinary = []
    for missing in combinations_with_replacement(range(len(tags)),deficit):
        if not any(tags[i] in guards for i in missing):
            continue
        levels = bounds.copy()
        for i in missing:
            levels[i] -= 1
        if min(levels) >= 0:
            ordinary.append(prod(shell(p,e) for p,e in zip(tags,levels)))
    weighted = []
    for levels in same_grade(sub_bounds,grades,rank):
        mass = guard_mass*prod(shell(p,e) for p,e in zip(raw+fine,levels))
        active_raw = prod(p**factors.get(p,0) for p,e in zip(raw,levels) if e)
        divisor = gcd(odd_part(mass),aux//active_raw)
        unit = odd_part(mass)//divisor
        weighted.append((unit,levels,mass))
    weighted.sort()
    return (tags,aux,K,background,sorted(ordinary,reverse=True),
            [mass for _,_,mass in weighted],[unit for unit,_,_ in weighted])


def main():
    args = partial_bounds()
    upper_norm, lower_norm = args[5:7]
    base = (1+upper_norm)*(3+lower_norm)/4
    for t in range(20):
        rank_bound = args[-2]+2*t
        ratio = rank_bound*(base/2+Fraction(3,2)*(1+upper_norm)*(1+lower_norm))/2**t
        assert ratio > Fraction(2,5)
    raw = [p for p in primerange(3,1000) if p not in args[0]][:7]
    safe = [p for p in primerange(3,100000) if p not in args[0]|set(raw)
            and set(factorint(odd_part(p-1))) <= set(raw)]
    data = spectator_data(args[0],raw,safe[:5],safe[5:11],[2]*5+[1]*2,2,8,3)
    result = verify_templates(args,*data[:-1],486,33,Fraction(5,2),
                              head_count=234,fine_units=data[-1])
    assert result == 133


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
