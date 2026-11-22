"""Exact checks for the 99-odd-prime three-stratum construction."""
import json
from fractions import Fraction
from itertools import product
from math import prod
from pathlib import Path
from sympy import isprime

import ppi_100_prime_pair as base


def main():
    folder = Path(__file__).parent
    tags = base.TAGS + (65537,)
    powers = base.POWERS + (1,)
    guards = [tags.index(p) for p in (*base.GUARDS, 65537)]
    orders = list(product(*(range(a+1) for a in powers)))
    masses = [base.shell_mass(tags, v) for v in orders]
    tag_order = prod(p**a for p, a in zip(tags, powers))
    assert sum(masses) == tag_order and isprime(65537)
    chosen = json.loads((folder / "ppi_99_slot_order.json").read_text())
    pos = {slot: j for j, slot in enumerate(chosen)}
    assert len(pos) == len(chosen) == 504
    guard = tuple(a if j in guards else 0 for j, a in enumerate(powers))
    guard_mass = base.shell_mass(tags, guard)
    background = {i for i,v in enumerate(orders) if i not in pos
                  and any(v[j] < powers[j] for j in guards)}
    stride = [prod(a+1 for a in powers[j+1:]) for j in range(len(tags))]
    for i, v in enumerate(orders):
        if v != guard and base.divides(v, guard):
            assert i in background
        if i not in pos and i not in background:
            continue
        if i in pos:
            assert any(v[j] < powers[j] for j in guards)
        for j, step in enumerate(stride):
            if v[j] == 0:
                continue
            parent = i-step
            if i in background:
                assert parent in background
            elif parent not in background:
                assert parent in pos and pos[parent] < pos[i]

    common = 83349
    variable = prod(p**a for p,a in zip(base.VARIABLE, base.VARIABLE_POWERS))
    control, fine = base.fine_slots(guard_mass, common)
    assert control == 34891
    M = min(v[3] & -v[3] for v in fine)
    assert M == 64*65536
    assert tag_order % common == 0
    K = tag_order//common * variable
    back = sum(masses[i] for i in background)
    assert back % common == 0 and all(masses[i] % common == 0 for i in chosen)
    n = back//common * variable
    weights = [masses[i]//common * variable for i in chosen]
    target = min(weights[:174])*99//100
    heights = [base.normalize(w,target) for w in weights]
    normalized = [h*w for h,w in zip(heights,weights)]
    groups = list(zip(normalized[::2],normalized[1::2]))
    fixed = base.bounds()
    axes,E,dy,upper,lower,U,L,eu,el,fixed_rank,_ = fixed
    assert fixed_rank == 46
    base.check_primes(axes)
    damping = base.PAD_PRIMES[:19]
    core = set(axes) | {q for p in damping for q in (p,p*(p-1)+1)}
    assert len(core | set(tags)) == 99 and not core.intersection(tags)
    rank = fixed_rank + 2*len(damping)
    margin = Fraction(27,10)
    slope = E*Fraction(round(Fraction(K*130352847*2**192,10**12*E)),2**192)
    delta = 2*slope/K
    depth = 0
    while 2*target > M << depth:
        depth += 1
    assert rank+depth+1 <= len(groups)
    groups = groups[:rank+depth+1]
    for i,h in enumerate(heights[:2*len(groups)]):
        assert 2*h <= 1 << max(0,i//2-rank)
    lows, highs = [], []
    frozen = 0
    for level in range(rank+1):
        minimum = frozen + sum(sum(g)/2**j for j,g in enumerate(groups[level:]))
        center = minimum + Fraction(5*M,4 << (rank-level)) - level*slope
        lows.append(center+slope/margin)
        highs.append(center+sum(groups[level])+target-slope/margin)
        frozen += 2*sum(groups[level])
    assert max(lows) < min(highs)
    intercept = (max(lows)+min(highs))/2
    c = 2*(n+intercept+rank*slope)/K
    assert 0 < c-rank*delta < c < 2 and 0 < delta < 1
    lower_helper = sorted(axes)[-4:]
    root = prod(q-1 for q in lower_helper)
    flat = root.bit_length()+min(lower_helper).bit_length()+1
    eta = Fraction(1,min(lower_helper)-1)
    alpha = Fraction((1 << flat)//root,1 << flat)
    assert 0 <= 1-alpha*root <= eta
    error = ((1+U)*(3+L)/4 + (1+L)*eta +
             Fraction(3,2)*rank*delta*(1+U)*(1+L))/(1 << len(damping))
    anchor_bound = 5+3*(1+upper)*(1+lower)
    blend_bits = 400
    blend = Fraction(((1 << blend_bits)//E)*E,1 << blend_bits)
    assert 0 < blend <= 1
    assert error + (1-blend)*anchor_bound < delta/margin
    def bits(value):
        d = value.denominator
        assert d & (d-1) == 0
        return d.bit_length()-1
    precision = (blend_bits+len(damping)+max(bits(K*c)+dy+flat+3,
                 max(bits(slope),bits(K*delta))+dy, eu+el)+6+M.bit_length())
    precision = max(precision,K.bit_length()+depth+rank+5+M.bit_length())
    assert precision <= len(fine)
    units = [v[0] for v in fine]
    extra = [(v[3] & -v[3]).bit_length()-M.bit_length() for v in fine]
    height = precision+max(extra)
    assert min(extra)==0 and units[:8]==[1]*8
    assert all(u > 0 and u%2 for u in units)
    frozen = 0
    for level in range(rank+1):
        offset = rank-level
        unit = Fraction(M,1 << offset)
        count = precision-offset
        assert count >= 8
        baseline = Fraction(sum(u << (count-j) for j,u in enumerate(units[:count],1)),1 << count)
        assert 0 <= baseline-Fraction(255,256) < Fraction(1,4)
        ordinary = [v/(1 << j) for j,g in enumerate(groups[level:]) for v in g]
        for i,h in enumerate(heights[:2*len(groups)]):
            assert bits(h/(1 << max(0,i//2-level))) <= height
        minimum = frozen + sum(ordinary) + unit*baseline
        for sign in (-1,0,1):
            wanted = intercept+level*slope+sign*slope/margin
            remaining = wanted-minimum-unit/4
            assert 0 <= remaining <= sum(groups[level])+target
            paid = minimum
            for value in ordinary:
                if value <= remaining:
                    remaining -= value
                    paid += value
            word = (remaining+unit/4)*(1 << count)/unit
            # Endpoints need not be grid points when margin has odd numerator.
            assert 1 << (count-2) <= word <= 3 << (count-2)
            if word.denominator != 1:
                assert sign != 0
                continue
            word = int(word)
            selected = 0
            for j in range(count,8,-1):
                if word%2:
                    word -= units[j-1]
                    selected += units[j-1] << (count-j)
                word //= 2
            assert 0 <= word < 256
            selected += word << (count-8)
            assert paid + unit*Fraction(selected,1 << count) == wanted
        frozen += 2*sum(groups[level])
    print(f"99 odd primes; ranks 0..{rank}; ordinary depth {depth}; "
          f"fine precision {precision}; dyadic exponent {height}.")
    print(f"Error including blend < 10/27 steps; {len(background)} background tags.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    main()

