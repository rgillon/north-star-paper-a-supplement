#!/usr/bin/env python3
"""Independent exact-arithmetic audit of PPI-H100, using data but no verifier imports.

The main certificate must also be run unchanged. This script independently checks
all prime certificates, donor incidence solutions and packet norms, the modulus,
spectator arithmetic, and the published precision/error budget. It is not a
replacement for the mathematical proof of annihilation or divisor closure.
"""
import ast
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
import sys
import time

if len(sys.argv) > 1:
    ROOT = Path(sys.argv[1]).resolve()
else:
    here = Path(__file__).resolve().parent
    candidates = [Path.cwd(), here, *here.parents]
    ROOT = next((p for p in candidates if (p / 'experiments' / 'ppi_100_certificate.json').is_file()
                 or (p / 'ppi_100_certificate.json').is_file()), None)
    if ROOT is None:
        raise SystemExit('Usage: python certificate_independent_check.py [repository_root]')
E = ROOT / 'experiments' if (ROOT / 'experiments' / 'ppi_100_certificate.json').is_file() else ROOT


def data(path):
    result = {}
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                result[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                pass
    return result


def v2(n):
    assert n > 0
    return (n & -n).bit_length() - 1


def phi(p, a):
    return 1 if a == 0 else (p - 1) * p ** (a - 1)


def mass(primes, levels):
    return math.prod(phi(p, a) for p, a in zip(primes, levels))


def solve(matrix, target):
    """Gauss-Jordan with exact fractions, independent of the SymPy inversion."""
    n = len(matrix[0])
    aug = [[F(x) for x in row] + [F(y)] for row, y in zip(matrix[:n], target[:n])]
    determinant = F(1)
    for col in range(n):
        pivot = next(row for row in range(col, n) if aug[row][col])
        if pivot != col:
            aug[pivot], aug[col] = aug[col], aug[pivot]
            determinant *= -1
        d = aug[col][col]
        determinant *= d
        aug[col] = [x / d for x in aug[col]]
        for row in range(n):
            if row != col:
                d = aug[row][col]
                aug[row] = [x - d * y for x, y in zip(aug[row], aug[col])]
    answer = [row[-1] for row in aug]
    assert all(sum(F(a) * b for a, b in zip(row, answer)) == rhs for row, rhs in zip(matrix, target))
    assert determinant.denominator == 1
    return answer, int(determinant)


def packet_norm(primes, masks, coefficients, denominator, lower=False):
    """Evaluate the defining packet sum directly at every activity pattern."""
    scaled = [int(denominator * c) for c in coefficients]
    assert all(denominator * c == a for c, a in zip(coefficients, scaled))
    rows = [[i for i in range(len(primes)) if m >> i & 1] for m in masks]
    if lower:
        reduced = []
        for a, row in zip(scaled, rows):
            q, r = divmod(a, math.prod(primes[i] - 1 for i in row))
            assert r == 0
            reduced.append(q)
    minimum = maximum = None
    for state in range(1 << len(primes)):
        value = -sum(scaled)
        for k, row in enumerate(rows):
            packet = 1
            for i in row:
                if lower:
                    packet *= -1 if state >> i & 1 else primes[i] - 1
                else:
                    packet *= 1 - primes[i] if state >> i & 1 else 1
            value += (reduced[k] if lower else scaled[k]) * packet
        minimum = value if minimum is None else min(minimum, value)
        maximum = value if maximum is None else max(maximum, value)
    return F(max(abs(minimum), abs(maximum)), denominator)


def sieve(bound):
    flags = bytearray([1]) * (bound + 1)
    flags[:2] = b'\0\0'
    for p in range(2, math.isqrt(bound) + 1):
        if flags[p]:
            flags[p*p::p] = b'\0' * ((bound - p*p) // p + 1)
    return [p for p, flag in enumerate(flags) if flag]


def pocklington(n, proven_factors):
    """Return a certificate using only already-proven primes dividing n-1."""
    known_product = 1
    witnesses = []
    for q in proven_factors:
        remainder, exponent = n - 1, 0
        while remainder % q == 0:
            remainder //= q
            exponent += 1
        assert exponent > 0
        known_product *= q ** exponent
        witness = next((a for a in range(2, 256)
                        if pow(a, n - 1, n) == 1
                        and math.gcd(pow(a, (n - 1) // q, n) - 1, n) == 1), None)
        assert witness is not None
        witnesses.append([str(q), exponent, witness])
    assert known_product ** 2 > n
    return {'prime': str(n), 'known_product': str(known_product), 'witnesses': witnesses}


def main():
    start = time.monotonic()
    h = data(E / 'ppi_100_prime_pair.py')
    u = data(E / 'ppi_unbalanced_pair.py')
    d = data(E / 'ppi_weighted_fibers.py')
    manifest = json.loads((E / 'ppi_100_certificate.json').read_text())
    powers = {int(p): a for p, a in manifest['odd_prime_powers'].items()}
    assert len(powers) == 100 and all(p > 2 and p % 2 and a > 0 for p, a in powers.items())
    odd_modulus = math.prod(p ** a for p, a in powers.items())
    modulus = odd_modulus << manifest['two_exponent']
    assert manifest['two_exponent'] == 3446
    assert modulus.bit_length() == 8123 and len(str(modulus)) == 2446

    # All small primes are certified by trial division, not a probable-prime test.
    small_primes = [p for p in powers if p < 1 << 64]
    divisors = sieve(math.isqrt(max(small_primes)))
    for p in small_primes:
        for q in divisors:
            if q * q > p:
                break
            assert p % q != 0
    large_primes = sorted(p for p in powers if p >= 1 << 64)
    assert len(small_primes) == 92 and len(large_primes) == 8
    certs = []
    for p in large_primes[:4]:
        remainder, factors = p - 1, []
        for q in divisors:
            if remainder % q == 0:
                factors.append(q)
                while remainder % q == 0:
                    remainder //= q
            if remainder == 1:
                break
        assert remainder == 1
        certs.append(pocklington(p, factors))
    for p in large_primes[4:]:
        certs.append(pocklington(p, [large_primes[0]]))

    # Reconstruct both donor interpolation systems directly from the data.
    hidden, hp, fibers, masks = u['HIDDEN'], u['HIDDEN_POWERS'], u['FIBERS'], u['MASKS']
    ps = [abs(row[0]) for row in fibers]
    states = list(itertools.product(*(range(a + 1) for a in reversed(hp))))
    states = [tuple(reversed(v)) for v in states]
    for p, c, exponent, bits in fibers:
        selected = [v for j, v in enumerate(states) if bits >> j & 1]
        assert sum(mass(hidden, v) for v in selected) == 2400 * c * (abs(p) - 1 if exponent == 1 else abs(p))
        for v in selected:
            assert all(w in selected for w in states if all(a <= b for a, b in zip(v, w)))
    assert sum((1 if p > 0 else -1) * c for p, c, _, _ in fibers) == 0
    assert all(sum(p * c for j, (p, c, _, _) in enumerate(fibers) if m >> j & 1) == 0 for m in masks)
    matrix = [[(m >> i) & 1 for m in masks] for i in range(len(ps))]
    coefficients, determinant = solve(matrix, [F(-1, p) for p in ps])
    assert determinant == 193
    du = determinant * math.prod(ps)
    U = packet_norm(ps, masks, coefficients, du)
    upper_initial = U
    assert U < 602953959432

    signed, lower_masks = d['LOWER']
    lp = list(map(abs, signed))
    assert sum((p - 1) * (1 if q > 0 else -1) for p, q in zip(lp, signed)) == 0
    assert all(sum(p for j, p in enumerate(signed) if m >> j & 1) == 0 for m in lower_masks)
    lm = [[(m >> i) & 1 for m in lower_masks] for i in range(len(lp))]
    lc, lower_determinant = solve(lm, [F(1 - p, p) for p in lp])
    dl = math.lcm(sum(lc).denominator, *(c.denominator * math.prod(p - 1 for i, p in enumerate(lp) if mask >> i & 1) for c, mask in zip(lc, lower_masks)))
    L = packet_norm(lp, lower_masks, lc, dl, lower=True)
    assert L < F(2447, 1000)
    odd = lambda n: n >> v2(n)
    Eodd = odd(du * dl)

    # Rebuild the eight large helpers; certify the modular interpolation residues.
    k = (249391 << 32) * math.prod(p for p in divisors if p <= 31) * du
    upper_helpers = [1 + k * a for a in (2, 12, 3, 4)]
    ol = odd(dl)
    P = upper_helpers[0]
    lower_helpers = [1 + P * (pow(P, -1, ol) + j * ol) for j in (304, 328, 170, 462)]
    assert set(upper_helpers + lower_helpers) == set(large_primes)
    assert sum(F(sign, p - 1) for sign, p in zip((1, 1, -1, -1), upper_helpers)) == 0
    assert sum(sign * (p - 1) for sign, p in zip((1, 1, -1, -1), lower_helpers)) == 0
    paired = []
    for qs, denom, upper in ((upper_helpers, du, True), (lower_helpers, ol, False)):
        norm_v = norm_c = 0
        for state in range(16):
            packet = math.prod(1 if state >> i & 1 else 1 - q for i, q in enumerate(qs))
            val = 1 - packet if upper else pow(-2, -1, ol) * (packet - math.prod(1 - q for q in qs))
            forced = state.bit_count() >= 3 if upper else state.bit_count() <= 1
            target = int(state.bit_count() == (3 if upper else 1))
            norm_v = max(norm_v, abs(val))
            if forced:
                assert (target - val) % denom == 0
                norm_c = max(norm_c, abs((target - val) // denom))
        norm = (du * U if upper else ol * L) * norm_c + norm_v
        denominator = du if upper else dl
        bits = (odd(denominator) * (int(norm) + 2)).bit_length()
        weight = F(((1 << bits) // odd(denominator)) * odd(denominator), 1 << bits)
        assert 0 < weight <= 1 and (1 - weight) * norm < 1
        paired.append(((U if upper else L) + 1, bits + v2(denominator)))
    dy = v2(du) + v2(dl)
    for primes in u['PAIRS']:
        Z = math.prod(p - 1 for p in primes)
        low, high = primes[0] - 1, primes[-1] - 1
        scale = min(F(2 * low * high, low + high), U * high / (U + high))
        bits = Z.bit_length() + 24
        alpha = F(((1 << bits) * scale) // Z, 1 << bits)
        error = max(abs(int(state.bit_count() == 1) + alpha * math.prod(1 if state >> i & 1 else 1 - q for i, q in enumerate(primes))) for state in range(1, 1 << len(primes)))
        U = max(F(1), alpha * Z, U * error)
        dy += bits
    eu, el = paired[0][1] + len(u['PAIRS']), paired[1][1]
    assert 30 < U < 31

    fixed = set(hidden + ps + lp + upper_helpers + lower_helpers + sum(u['PAIRS'], []))
    damping = {p for q in h['PAD_PRIMES'] for p in (q, q * (q - 1) + 1)}
    tags = h['TAGS']
    assert len(fixed) == 48 and len(damping) == 40 and len(tags) == 12
    assert len(fixed | damping | set(tags)) == 100
    reconstructed = dict.fromkeys(fixed, 1)
    reconstructed.update(zip(hidden, hp))
    reconstructed.update((abs(p), a) for p, _, a, _ in fibers)
    for p in h['PAD_PRIMES']:
        reconstructed.update({p: 2, p * (p - 1) + 1: 1})
    reconstructed.update(zip(tags, h['POWERS']))
    for p, a in zip(h['VARIABLE'], h['VARIABLE_POWERS']):
        reconstructed[p] += a
    for p in h['CONTROL']:
        reconstructed[p] += 1
    assert reconstructed == powers

    # Independent spectator mass data and fine-weight arithmetic.
    guards = h['GUARDS']
    guard_tuple = tuple(a if p in guards else 0 for p, a in zip(tags, h['POWERS']))
    G = mass(tags, guard_tuple)
    C = math.prod(h['CONTROL'])
    common = 83349
    assert G == 64 * common and C == 34891
    Q = math.prod(p ** a for p, a in zip(tags, h['POWERS']))
    B = math.prod(p ** a for p, a in zip(h['VARIABLE'], h['VARIABLE_POWERS']))
    K = Q * B // common
    fine = []
    for tup in itertools.product(*(range(a + 1) for a in h['VARIABLE_POWERS'])):
        V = G * mass(h['VARIABLE'], tup)
        index = math.prod(p for p in h['CONTROL'] if V % p == 0)
        assert math.gcd(V, C * C) == index
        assert (V >> v2(V)) % (common * index) == 0
        unit = (V >> v2(V)) // (common * index)
        fine.append((unit, sum(tup), tup, V, index))
    fine.sort()
    assert len(fine) == 3840 and all(row[0] == 1 for row in fine[:8])
    assert fine[8][0] != 1
    extra = max(v2(row[3]) for row in fine) - 6
    assert extra == 7

    tuples = list(itertools.product(*(range(a + 1) for a in h['POWERS'])))
    masses = [mass(tags, t) for t in tuples]
    chosen = json.loads((E / 'ppi_100_slot_order.json').read_text())
    early, late, background = [], [], []
    ordinary = set(chosen)
    gi = [tags.index(p) for p in guards]
    for i, t in enumerate(tuples):
        missing = [t[j] < h['POWERS'][j] for j in gi]
        stratum = 0 if any(missing[:2]) else (1 if missing[2] else 2)
        if stratum < 2:
            if i not in ordinary:
                background.append(i)
            if masses[i] % common == 0:
                (early if stratum == 0 else late).append(i)
    early.sort(key=lambda i: (masses[i], tuples[i]), reverse=True)
    late.sort(key=lambda i: (masses[i], tuples[i]), reverse=True)
    assert set(chosen[:214]) == set(early[:214]) and set(chosen[214:]) == set(late[:400])
    assert len(chosen) == len(ordinary) == 614
    n = sum(masses[i] for i in background) * B // common
    assert sum(masses[i] for i in background) % common == 0
    target = (masses[early[213]] * B // common) * 99 // 100
    weights = [masses[i] * B // common for i in chosen]
    heights = []
    for weight in weights:
        height = F(1)
        while weight * height > target:
            height /= 2
        while 2 * weight * height <= target:
            height *= 2
        heights.append(height)
    normalized = [w * a for w, a in zip(weights, heights)]
    groups = [sum(normalized[i:i+2]) for i in range(0, 614, 2)]
    step = F(int(manifest['step_numerator']), 1 << manifest['step_binary_exponent'])
    intercept = F(int(manifest['intercept_numerator']), 1 << manifest['intercept_binary_exponent'])
    s = (Eodd * (1 << 20) // target).bit_length()
    assert s == 161
    assert step == F(Eodd * round(F(2 * sum(groups[:107]), 107) * (1 << s) / Eodd), 1 << s)
    R, depth, M, m = 86, 159, 64, F(7, 2)
    delta = 2 * step / K
    c = 2 * (n + intercept + R * step) / K
    assert 0 < c - R * delta < c < 2
    eta = F(1, min(lower_helpers) - 1)
    rational_ratio = (((1 + U) * (3 + L) / 4 + (1 + L) * eta) / delta + F(3, 2) * R * (1 + U) * (1 + L)) / (1 << 20)
    assert rational_ratio < F(254, 1000)
    anchor_difference = 5 + 3 * (1 + paired[0][0]) * (1 + paired[1][0])
    blend_bits = (33 * m * Eodd * anchor_difference // delta).bit_length() + 1
    assert blend_bits == 367
    lam = F(((1 << blend_bits) // Eodd) * Eodd, 1 << blend_bits)
    blend_ratio = (1 - lam) * anchor_difference / delta
    assert blend_ratio < 1 / (33 * m)
    assert rational_ratio + blend_ratio < 1 / m
    flat_bits = math.prod(p - 1 for p in lower_helpers).bit_length() + min(lower_helpers).bit_length() + 1
    base_bits = v2((K * c).denominator)
    assert (K * c).denominator == 1 << base_bits
    assert (flat_bits, base_bits) == (2003, 177)
    profile_bits = blend_bits + 20 + max(base_bits + dy + flat_bits + 3, s + dy, eu + el) + 6 + M.bit_length()
    ordinary_bits = K.bit_length() + depth + R + 5 + M.bit_length()
    A = max(profile_bits, ordinary_bits)
    assert A == 3353 and A + extra + R == 3446
    tail = sum(F(row[0], 1 << j) for j, row in enumerate(fine[8:A], 9))
    assert tail < F(1, 4)
    lows, highs = [], []
    frozen = F(0)
    for level in range(R + 1):
        C_l = frozen + sum(F(g, 1 << j) for j, g in enumerate(groups[level:R+depth+1]))
        nominal = C_l + F(5 * M, 4 << (R - level)) - level * step
        lows.append(nominal + step / m)
        highs.append(nominal + groups[level] + target - step / m)
        frozen += 2 * groups[level]
    assert max(lows) < min(highs)
    assert intercept == (max(lows) + min(highs)) / 2

    summary = {
        'scope': 'Independent arithmetic cross-check; symbolic proof reviewed separately.',
        'source_commit': '09d77232d96df1e16ef52917ddf46ccff90deb0b',
        'result': 'PASS',
        'python': sys.version,
        'seconds': time.monotonic() - start,
        'prime_certificates': {'trial_division_count': len(small_primes), 'pocklington_count': len(certs), 'largest_trial_divisor': divisors[-1]},
        'modulus': {'odd_distinct': 100, 'all_distinct': 101, 'odd_Omega': sum(powers.values()), 'bit_length': 8123, 'decimal_digits': 2446, 'two_exponent': 3446},
        'donor': {'upper_minor': determinant, 'lower_minor': lower_determinant, 'upper_denominator': str(du), 'lower_denominator': str(dl), 'E': str(Eodd), 'dy': dy, 'eu': eu, 'el': el, 'upper_initial_norm': float(upper_initial), 'upper_norm': float(U), 'lower_norm': float(L)},
        'spectator': {'tags': len(tuples), 'background': len(background), 'ordinary_slots': len(chosen), 'fine_slots': len(fine), 'fine_used': A, 'levels': R + 1, 'K': str(K), 'n': str(n), 'H0': str(target), 'control_kernel': C, 'common': common, 'fine_tail': float(tail)},
        'precision': {'base': base_bits, 'step_grid': s, 'flat': flat_bits, 'blend': blend_bits, 'profile': profile_bits, 'ordinary': ordinary_bits, 'A': A, 'extra': extra, 'R': R, 'dyadic': A + extra + R},
        'errors': {'rational_ratio': float(rational_ratio), 'blend_ratio': float(blend_ratio), 'total_ratio': float(rational_ratio + blend_ratio), 'allowed_ratio': float(1 / m), 'interval_width_relative_step': float((min(highs) - max(lows)) / step)},
        'exact_rationals': {name: {'numerator': str(value.numerator), 'denominator': str(value.denominator)} for name, value in {
            'upper_norm': U, 'lower_norm': L, 'rational_error_over_delta': rational_ratio,
            'blend_error_over_delta': blend_ratio, 'total_error_over_delta': rational_ratio + blend_ratio,
            'allowed_error_over_delta': 1 / m, 'fine_tail': tail,
            'interval_width_relative_step': (min(highs) - max(lows)) / step,
            'step': step, 'intercept': intercept,
        }.items()},
        'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (E / 'ppi_100_prime_pair.py', E / 'ppi_100_certificate.json', E / 'ppi_100_slot_order.json', E / 'ppi_unbalanced_pair.py', E / 'ppi_weighted_fibers.py')},
        'pocklington': certs,
    }
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('Run with assertions enabled.')
    main()
