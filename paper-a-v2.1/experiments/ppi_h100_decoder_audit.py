"""Independent standard-library spectator/decoder audit for PPI-H100.

No imports from the construction or its verification modules.  The two JSON
files at the pinned source commit are the only construction input files.
Checks all tag cover edges, all fine-slot cover edges, all 87 rank levels,
the common intercept, all physical subgroup-height inequalities, and a
deterministic selection of interior mass-grid values.  The accompanying
audit note proves coverage of every mass-grid value, not only these samples.
"""
from fractions import Fraction as F
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import random
import sys

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
# Accept either a repository root or a flat certificate-bundle directory.
DATA = ROOT / 'experiments' if (ROOT / 'experiments/ppi_100_certificate.json').is_file() else ROOT
TAGS = (3, 7, 19, 23, 37, 41, 43, 47, 71, 83, 113, 149)
TAG_EXP = (5, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1)
GUARDS = (3, 7, 113)
CONTROL = (23, 37, 41)
VARIABLE = (19, 43, 47, 71, 83, 149)
VAR_EXP = (4, 3, 3, 3, 3, 2)
R, A, DEPTH, DYADIC, M = 86, 3353, 159, 3446, 64

def phi(p, a):
    return 1 if a == 0 else (p - 1) * p ** (a - 1)

def mass(ps, es):
    return prod(phi(p, e) for p, e in zip(ps, es))

def valuation(n, p):
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e

def parent_orders(v):
    for j, e in enumerate(v):
        if e:
            yield v[:j] + (e - 1,) + v[j + 1:]

def dyadic(n, exponent):
    return F(n, 1 << exponent) if exponent >= 0 else F(n << -exponent)

def integer(x):
    assert x.denominator == 1
    return x.numerator

def main():
    mf_path = DATA / 'ppi_100_certificate.json'
    slots_path = DATA / 'ppi_100_slot_order.json'
    manifest = json.loads(mf_path.read_text())
    chosen = json.loads(slots_path.read_text())
    step = dyadic(int(manifest['step_numerator']), manifest['step_binary_exponent'])
    intercept = dyadic(int(manifest['intercept_numerator']), manifest['intercept_binary_exponent'])
    assert manifest['two_exponent'] == DYADIC
    orders = list(product(*(range(e + 1) for e in TAG_EXP)))
    idx = {v: j for j, v in enumerate(orders)}
    weights = [mass(TAGS, v) for v in orders]
    Q = prod(p ** e for p, e in zip(TAGS, TAG_EXP))
    B = prod(p ** e for p, e in zip(VARIABLE, VAR_EXP))
    C = prod(CONTROL)
    h = 3 ** 5 * 7 ** 3
    guard = tuple(e if p in GUARDS else 0 for p, e in zip(TAGS, TAG_EXP))
    G = mass(TAGS, guard)
    assert (len(orders), sum(weights), G, C) == (124416, Q, 64 * h, 34891)
    K = Q * B // h
    assert K * h == Q * B
    for p, a in zip(TAGS, TAG_EXP):
        extra = VAR_EXP[VARIABLE.index(p)] if p in VARIABLE else (1 if p in CONTROL else 0)
        assert manifest['odd_prime_powers'][str(p)] == a + extra

    early_all = [j for j, v in enumerate(orders) if v[0] < 5 or v[1] < 3]
    late_all = [j for j, v in enumerate(orders) if v[0] == 5 and v[1] == 3 and v[10] == 0]
    early = sorted((j for j in early_all if weights[j] % h == 0), key=lambda j: (weights[j], j), reverse=True)[:214]
    late = sorted((j for j in late_all if weights[j] % h == 0), key=lambda j: (weights[j], j), reverse=True)[:400]
    assert len(chosen) == len(set(chosen)) == 614
    assert set(chosen[:214]) == set(early) and set(chosen[214:]) == set(late)
    pos = {j: k for k, j in enumerate(chosen)}
    background = set(early_all + late_all) - set(chosen)
    assert len(background) == 121210
    background_sum = sum(weights[j] for j in background)
    assert background_sum % h == 0
    n = background_sum * B // h

    tag_edges = 0
    ordinary_edges = []
    for j, v in enumerate(orders):
        if v != guard and all(a <= b for a, b in zip(v, guard)):
            assert j in background
        for parent in parent_orders(v):
            q = idx[parent]
            tag_edges += 1
            if j in background:
                assert q in background
            if j in pos:
                assert q in background or (q in pos and pos[q] < pos[j])
                if q in pos:
                    assert weights[j] >= 2 * weights[q]
                    ordinary_edges.append((pos[q], pos[j]))

    W = [weights[j] * B // h for j in chosen]
    H = 99 * min(W[:214]) // 100
    beta_exp = []
    z = []
    for w in W:
        exp = 0
        val = F(w)
        while val > H:
            exp -= 1
            val /= 2
        while val * 2 <= H:
            exp += 1
            val *= 2
        assert H < 2 * val <= 2 * H
        assert val == w * F(2) ** exp
        beta_exp.append(exp)
        z.append(val)
    for parent, child in ordinary_edges:
        assert beta_exp[child] <= beta_exp[parent] - 1
    groups = list(zip(z[::2], z[1::2]))
    used = groups[:R + DEPTH + 1]
    assert len(used) == 246 and H <= (M << DEPTH) // 2

    fine = []
    for v in product(*(range(e + 1) for e in VAR_EXP)):
        raw = G * mass(VARIABLE, v)
        two = valuation(raw, 2)
        control_index = prod(p ** valuation(raw, p) for p in CONTROL)
        assert C % control_index == 0
        unit, rem = divmod(raw >> two, h * control_index)
        assert rem == 0 and unit > 0 and unit % 2 == 1
        fine.append((unit, sum(v), v, raw, control_index, two))
    fine.sort()
    fpos = {x[2]: j for j, x in enumerate(fine)}
    fine_edges = []
    for j, f in enumerate(fine):
        for pv in parent_orders(f[2]):
            i = fpos[pv]
            p = fine[i]
            assert i < j and f[3] % p[3] == 0 and f[4] % p[4] == 0
            assert f[0] % p[0] == 0 and p[5] <= f[5]
            if j < A:
                fine_edges.append((i, j))
    assert len(fine) == 3840
    assert [f[0] for f in fine[:8]] == [1] * 8 and fine[8][0] > 1
    tail = sum(fine[j - 1][0] << (A - j) for j in range(9, A + 1))
    assert 4 * tail < 1 << A
    fine_total = sum(fine[j - 1][0] << (A - j) for j in range(1, A + 1))
    T = F(fine_total, 1 << A)
    assert max(f[5] for f in fine) - 6 == 7

    def ordinary_height_range(slot, level):
        k = slot // 2
        if k >= len(used):
            return None
        if k < level:
            t = DYADIC + beta_exp[slot] + 1
            return (t, t)
        lo = DYADIC + beta_exp[slot] - (k - level)
        return (lo, lo + 1)

    def fine_height_range(slot, level):
        # slot is zero-based; the first nominal fine height has j = 1.
        lo = DYADIC + level - R - (slot + 1) + 6 - fine[slot][5]
        return (lo, lo + 1)

    ordinary_height_checks = fine_height_checks = 0
    for level in range(R + 1):
        for j in range(614):
            now = ordinary_height_range(j, level)
            if now:
                assert 0 <= now[0] <= now[1] <= DYADIC
                ordinary_height_checks += 1
                if level:
                    prev = ordinary_height_range(j, level - 1)
                    assert now[0] >= prev[1]
        for parent, child in ordinary_edges:
            ch = ordinary_height_range(child, level)
            if ch:
                pa = ordinary_height_range(parent, level)
                assert pa is not None and pa[0] >= ch[1]
        for j in range(A):
            now = fine_height_range(j, level)
            assert 0 <= now[0] <= now[1] <= DYADIC
            fine_height_checks += 1
            if level:
                assert now[0] >= fine_height_range(j, level - 1)[1]
        for parent, child in fine_edges:
            assert fine_height_range(parent, level)[0] >= fine_height_range(child, level)[1]

    lows, highs = [], []
    minimums = []
    frozen = F(0)
    for level in range(R + 1):
        ordinary_min = frozen + sum(sum(g) / (1 << j) for j, g in enumerate(used[level:]))
        minimums.append(ordinary_min)
        fine_scale = F(M, 1 << (R - level))
        nominal = ordinary_min + 5 * fine_scale / 4 - level * step
        lows.append(nominal + 2 * step / 7)
        highs.append(nominal + sum(used[level]) + H - 2 * step / 7)
        frozen += 2 * sum(used[level])
    intersection_lo, intersection_hi = max(lows), min(highs)
    assert intersection_lo < intersection_hi
    assert intercept == (intersection_lo + intersection_hi) / 2

    # Alternate implementation of the fine decoder: accumulate a triangular
    # subset sum in one integer, extracting its forced bits without repeated
    # subtraction-and-division of the working word used in the source code.
    def decode_fine(q):
        selected = 0
        for j in range(A, 8, -1):
            shift = A - j
            if ((q - selected) >> shift) & 1:
                selected += fine[j - 1][0] << shift
        high, residue = divmod(q - selected, 1 << (A - 8))
        assert residue == 0 and 0 <= high < 256
        reconstructed = selected + sum(1 << (A - j) for j in range(1, 9) if high & (1 << (8 - j)))
        assert reconstructed == q
        return reconstructed

    scale = 1 << DYADIC
    rng = random.Random(1003446)
    target_checks = 0
    smallest_endpoint_margin = None
    for level in range(R + 1):
        fine_scale = F(M, 1 << (R - level))
        unit = integer(fine_scale * scale)
        tick = unit >> A
        assert tick > 0 and unit == tick << A
        ordinary_min = integer(minimums[level] * scale)
        baseline = n * scale + ordinary_min + fine_total * tick
        reserve = unit // 4
        center = integer((n + intercept + level * step) * scale)
        radius = integer(2 * step * scale / 7)
        assert center % tick == radius % tick == 0
        low, high = center - radius, center + radius
        cap = integer((sum(used[level]) + H) * scale)
        for bound in (low, high):
            residual = bound - baseline - reserve
            assert 0 <= residual <= cap
            gap = F(min(residual, cap - residual), scale)
            smallest_endpoint_margin = gap if smallest_endpoint_margin is None else min(gap, smallest_endpoint_margin)
        candidates = {low, low + tick, center, high - tick, high}
        for _ in range(20):
            candidates.add(rng.randrange(low // tick, high // tick + 1) * tick)
        active = [integer(z * scale / (1 << j)) for j, pair in enumerate(used[level:]) for z in pair]
        for requested in candidates:
            residual = requested - baseline - reserve
            decoded = baseline
            for amount in active:
                if amount <= residual:
                    residual -= amount
                    decoded += amount
            q, rem = divmod(residual + reserve, tick)
            assert rem == 0 and (1 << (A - 2)) <= q <= 3 << (A - 2)
            decoded += decode_fine(q) * tick
            assert decoded == requested
            target_checks += 1

    result = {
        'status': 'PASS',
        'scope': 'Independent finite spectator/decoder premises; universal coverage proof in audit note',
        'source_sha256': {'experiments/' + p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (mf_path, slots_path)},
        'tag_orders': len(orders),
        'tag_cover_edges_checked': tag_edges,
        'background_tags': len(background),
        'ordinary_slots': len(chosen),
        'ordinary_cover_edges': len(ordinary_edges),
        'ordinary_groups_used': len(used),
        'fine_slots': len(fine),
        'fine_slots_used': A,
        'used_fine_cover_edges': len(fine_edges),
        'levels': R + 1,
        'ordinary_height_ranges_checked': ordinary_height_checks,
        'fine_height_ranges_checked': fine_height_checks,
        'independent_exact_targets_replayed': target_checks,
        'fine_prefix_unit_weights': 8,
        'fine_tail_decimal': float(F(tail, 1 << A)),
        'dyadic_exponent': DYADIC,
        'smallest_used_fine_subgroup_exponent': min(fine_height_range(j, 0)[0] for j in range(A)),
        'Q': str(Q), 'B': str(B), 'C': C, 'h': h, 'G': G,
        'K': str(K), 'n': str(n), 'ordinary_target': str(H),
        'intercept_intersection_width_over_step': str((intersection_hi - intersection_lo) / step),
        'minimum_endpoint_slack_over_step_decimal': float(smallest_endpoint_margin / step),
    }
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Assertions must be enabled')
    main()
