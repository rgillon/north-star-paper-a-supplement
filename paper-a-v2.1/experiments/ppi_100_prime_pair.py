"""Exact certificate for two nontranslated homometric subgroup unions on 100 odd primes."""
import json
from fractions import Fraction
from itertools import product
from math import isqrt, prod
from pathlib import Path

from ppi_mixed_slots import verify_templates
from ppi_rank_box_pair import odd_part
from ppi_unbalanced_pair import FIBERS, HIDDEN, HIDDEN_POWERS, bounds
from ppi_weighted_digits import shell


TAGS = (3, 7, 19, 23, 37, 41, 43, 47, 71, 83, 113, 149)
POWERS = (5, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1)
GUARDS = (3, 7, 113)
CONTROL = (23, 37, 41)
VARIABLE = (19, 43, 47, 71, 83, 149)
VARIABLE_POWERS = (4, 3, 3, 3, 3, 2)
PAD_PRIMES = (2097397, 2097481, 2097523, 2097673, 2097727, 2098081, 2098279,
              2098471, 2099017, 2099089, 2099479, 2099593, 2099659, 2100001,
              2100463, 2100841, 2101003, 2101243, 2101261, 2101327)
EARLY, LATE = 214, 400
MARGIN = Fraction(7, 2)


def divides(v, w):
    return all(a <= b for a, b in zip(v, w))


def shell_mass(primes, order):
    return prod(shell(p, a) for p, a in zip(primes, order))


def control_part(n):
    result = 1
    for p in CONTROL:
        while n % p == 0:
            n //= p
            result *= p
    return result


def check_primes(fixed):
    pads = {p for r in PAD_PRIMES for p in (r, r*(r-1)+1)}
    assert len(fixed) == 48 and len(pads) == 40 and len(set(TAGS)) == 12
    assert len(fixed | pads | set(TAGS)) == 100
    # Independent trial-division certificate for every listed damping and tag prime.
    limit = isqrt(max(pads))
    sieve = bytearray([1])*(limit+1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(limit)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit-p*p)//p+1)
    divisors = [p for p, prime in enumerate(sieve) if prime]
    for q in pads | set(TAGS):
        assert q > 2 and all(q % p for p in divisors if p*p <= q)


def fine_slots(guard_mass, common):
    """All variable kernel shells, their odd weights and control-subgroup indices."""
    kernel = control_part(guard_mass*shell_mass(VARIABLE, VARIABLE_POWERS))
    slots = []
    for v in product(*(range(a+1) for a in VARIABLE_POWERS)):
        mass = guard_mass*shell_mass(VARIABLE, v)
        index = control_part(mass)
        assert kernel % index == 0 and odd_part(mass) % (common*index) == 0
        unit = odd_part(mass)//(common*index)
        slots.append((unit, sum(v), v, mass, index))
    slots.sort()
    position = {slot[2]: j for j, slot in enumerate(slots)}
    for j, (_, _, v, mass, index) in enumerate(slots):
        for i, a in enumerate(v):
            if a:
                parent = position[v[:i]+(a-1,)+v[i+1:]]
                assert parent < j and index % slots[parent][4] == 0
                assert mass % slots[parent][3] == 0
    return kernel, slots


def normalize(weight, target):
    """The unique power-of-two height sending weight into (target/2, target]."""
    shift = weight.bit_length()-target.bit_length()
    nominal = Fraction(2)**(-shift)
    if nominal*weight > target:
        nominal /= 2
    assert target < 2*nominal*weight <= 2*target
    return nominal



def check_modulus(manifest, fixed, control_kernel):
    powers = dict.fromkeys(fixed, 1)
    powers.update(zip(HIDDEN, HIDDEN_POWERS))
    powers.update((abs(p), a) for p, _, a, _ in FIBERS)
    for p in PAD_PRIMES:
        powers.update({p: 2, p*(p-1)+1: 1})
    powers.update(zip(TAGS, POWERS))
    for p, a in zip(VARIABLE, VARIABLE_POWERS):
        powers[p] += a
    for p in CONTROL:
        while control_kernel % p == 0:
            powers[p] += 1
            control_kernel //= p
    assert control_kernel == 1
    assert powers == {int(p): a for p, a in manifest["odd_prime_powers"].items()}
    assert manifest["two_exponent"] == 3446


def replay_mass_grid(manifest, groups, fine, K, n, target, fixed):
    """Reconstruct physical masses as integers at full dyadic precision."""
    step = Fraction(int(manifest["step_numerator"]), 1 << manifest["step_binary_exponent"])
    intercept = Fraction(int(manifest["intercept_numerator"]), 1 << manifest["intercept_binary_exponent"])
    E = fixed[1]
    mean = sum(map(sum, groups[:EARLY//2]))/(EARLY//2)
    assert step == Fraction(E*round(2*mean*(1 << 161)/E), 1 << 161)
    scale = 1 << manifest["two_exponent"]
    depth, digits, rank = 159, 3353, 86
    M = min(x[3] & -x[3] for x in fine)
    assert M == 64 and M*(1 << (depth-1)) < 2*target <= M*(1 << depth)
    assert all(x[0] == 1 for x in fine[:8])
    assert 4*sum(x[0] << (digits-j) for j, x in enumerate(fine[8:digits], 9)) < 1 << digits
    delta = 2*step/K
    c = 2*(n+intercept+rank*step)/K
    U, L = fixed[5:7]
    eta = Fraction(1, min(sorted(fixed[0])[-4:])-1)
    base = (1+U)*(3+L)/4+(1+L)*eta
    ratio = (base/delta+Fraction(3,2)*rank*(1+U)*(1+L))/(1 << 20)
    assert 0 < c-rank*delta < c < 2 and ratio < Fraction(254,1000)
    assert ratio+1/(33*MARGIN) < 1/MARGIN
    # The stored intercept is the exact midpoint in the defining construction.
    lo, hi = [], []
    used = groups[:rank+depth+1]
    frozen = 0
    for level in range(rank+1):
        active = [v/(1 << j) for j, pair in enumerate(used[level:]) for v in pair]
        unit = Fraction(M, 1 << (rank-level))
        nominal = frozen+sum(active)+5*unit/4-level*step
        lo.append(nominal+step/MARGIN)
        hi.append(nominal+sum(used[level])+target-step/MARGIN)
        frozen += 2*sum(used[level])
    assert max(lo) < min(hi) and intercept == (max(lo)+min(hi))/2

    def integral(value):
        assert value.denominator == 1
        return int(value)

    frozen = 0
    for level in range(rank+1):
        tick = (M*scale) >> (rank-level+digits)
        active = [integral(v*scale/(1 << j))
                  for j, pair in enumerate(used[level:]) for v in pair]
        fine_unit = M*scale >> (rank-level)
        fine_base = sum(x[0]*(fine_unit >> j) for j, x in enumerate(fine[:digits], 1))
        minimum = n*scale+frozen+sum(active)+fine_base
        reserve = fine_unit//4
        center = integral((n+intercept+level*step)*scale)
        radius = integral(step*scale/MARGIN)
        for offset in (-radius, -radius+tick, 0, radius-tick, radius):
            requested = center+offset
            left = requested-minimum-reserve
            assert 0 <= left <= integral((sum(used[level])+target)*scale)
            actual = minimum
            for amount in active:
                if left >= amount:
                    left -= amount
                    actual += amount
            word, rem = divmod(left+reserve, tick)
            assert rem == 0 and (1 << (digits-2)) <= word <= 3 << (digits-2)
            for j in range(digits, 8, -1):
                weight = fine[j-1][0]
                if word & 1:
                    word -= weight
                    actual += weight*(fine_unit >> j)
                word //= 2
            assert 0 <= word < 256
            for j in range(1, 9):
                if word & (1 << (8-j)):
                    actual += fine_unit >> j
            assert actual == requested
        frozen += integral(2*sum(used[level])*scale)
    print("Integer mass replay passed 435 targets across 87 levels. Rational error < 0.254 steps.")


def main():
    fixed = bounds()
    check_primes(fixed[0])
    assert set(TAGS) == set(GUARDS+CONTROL+VARIABLE)
    guard_ids = tuple(TAGS.index(p) for p in GUARDS)
    guard_order = tuple(a if i in guard_ids else 0 for i, a in enumerate(POWERS))
    guard_mass = shell_mass(TAGS, guard_order)
    common = odd_part(guard_mass)//control_part(guard_mass)
    quotient = prod(p**a for p, a in zip(TAGS, POWERS))
    variable_kernel = prod(p**a for p, a in zip(VARIABLE, VARIABLE_POWERS))
    assert common == 83349 and quotient % common == 0
    K = quotient*variable_kernel//common
    control_kernel, fine = fine_slots(guard_mass, common)

    orders = list(product(*(range(a+1) for a in POWERS)))
    masses = [shell_mass(TAGS, v) for v in orders]
    assert sum(masses) == quotient
    strata = []
    early, late = [], []
    for i, (v, mass) in enumerate(zip(orders, masses)):
        missing = [v[j] < POWERS[j] for j in guard_ids]
        stratum = 0 if any(missing[:2]) else 1 if missing[2] else 2
        strata.append(stratum)
        if stratum != 2 and mass % common == 0:
            (early if stratum == 0 else late).append(i)
    # For equal weights the larger lexicographic order comes first.
    early = sorted(early, key=lambda i: (masses[i], i), reverse=True)[:EARLY]
    late = sorted(late, key=lambda i: (masses[i], i), reverse=True)[:LATE]
    chosen = json.loads(Path(__file__).with_name("ppi_100_slot_order.json").read_text())
    assert len(chosen) == len(set(chosen)) == EARLY+LATE
    assert set(chosen[:EARLY]) == set(early) and set(chosen[EARLY:]) == set(late)
    position = {index: j for j, index in enumerate(chosen)}
    # Check background closure, every strict subtag of the fine tag, and the
    # direction of every cover edge between ordinary slots independently.
    background = {i for i, k in enumerate(strata) if k != 2 and i not in position}
    stride = [prod(a+1 for a in POWERS[j+1:]) for j in range(len(TAGS))]
    for i, v in enumerate(orders):
        if v != guard_order and divides(v, guard_order):
            assert i in background
        for j, a in enumerate(v):
            if not a:
                continue
            parent = i-stride[j]
            if i in background:
                assert parent in background
            if i in position and parent not in background:
                assert parent in position and position[parent] < position[i]
                assert masses[i] >= 2*masses[parent]
    background_mass = sum(masses[i] for i in background)
    assert background_mass % common == 0
    n = background_mass*variable_kernel//common
    ordinary = [masses[i]*variable_kernel//common for i in early+late]
    target = ordinary[EARLY-1]*99//100
    ordered_weights = [masses[i]*variable_kernel//common for i in chosen]
    heights = [normalize(w, target) for w in ordered_weights]
    assert all(2*h <= 1 for h in heights[:EARLY])
    assert all(2*h <= 2**(EARLY//2-86) for h in heights[EARLY:])
    normalized = [h*w for h, w in zip(heights, ordered_weights)]
    groups = list(zip(normalized[::2], normalized[1::2]))
    manifest = json.loads(Path(__file__).with_name("ppi_100_certificate.json").read_text())
    check_modulus(manifest, fixed[0], control_kernel)
    replay_mass_grid(manifest, groups, fine, K, n, target, fixed)
    total = verify_templates(fixed, list(TAGS), control_kernel, K, n, ordinary,
                             [x[3] for x in fine], EARLY+LATE, len(PAD_PRIMES), MARGIN,
                             head_count=EARLY, fine_units=[x[0] for x in fine],
                             prepared_groups=(target, groups))
    assert total == 100
    assert len(orders) == 124416 and len(fine) == 3840
    print(f"Tag orders={len(orders)}, full background tags={len(background)}, "
          f"common odd divisor={common}, control kernel={control_kernel}.")
    print(f"Normalized total K={K}, background n={n}, ordinary target={target}.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
