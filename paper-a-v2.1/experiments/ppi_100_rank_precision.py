"""Verify the rank-dependent truncation of the H100 fine-slot decoder."""
from fractions import Fraction
from itertools import product
from math import prod
import json
from pathlib import Path

from ppi_100_prime_pair import (
    EARLY, MARGIN, POWERS, TAGS, VARIABLE, VARIABLE_POWERS,
    fine_slots, normalize, shell_mass,
)


def main():
    folder = Path(__file__).parent
    manifest = json.loads((folder / "ppi_100_certificate.json").read_text())
    chosen = json.loads((folder / "ppi_100_slot_order.json").read_text())
    orders = list(product(*(range(a + 1) for a in POWERS)))
    variable = prod(p**a for p, a in zip(VARIABLE, VARIABLE_POWERS))
    common = 83349
    weights = []
    for pos in chosen:
        mass = shell_mass(TAGS, orders[pos])
        assert mass % common == 0
        weights.append(mass // common * variable)
    target = min(weights[:EARLY]) * 99 // 100
    heights = [normalize(w, target) for w in weights]
    normalized = [h * w for h, w in zip(heights, weights)]
    groups = list(zip(normalized[::2], normalized[1::2]))[:246]
    assert len(groups) == 246
    step = Fraction(int(manifest["step_numerator"]), 2**manifest["step_binary_exponent"])
    intercept = Fraction(int(manifest["intercept_numerator"]), 2**manifest["intercept_binary_exponent"])
    _, slots = fine_slots(5334336, common)
    units = [s[0] for s in slots]
    extra = [(s[3] & -s[3]).bit_length() - 7 for s in slots]
    assert min(extra) == 0 and max(extra) == 7 and units[:8] == [1] * 8
    rank, precision, height = 86, 3353, 3360
    max_used = 0
    frozen = 0
    for level in range(rank + 1):
        offset = rank - level
        digits = precision - offset
        unit = Fraction(64, 1 << offset)
        for j in range(1, digits + 1):
            max_used = max(max_used, offset + j + extra[j-1])
        fine_integer = sum(w << (digits-j) for j, w in enumerate(units[:digits], 1))
        baseline = Fraction(fine_integer, 1 << digits)
        assert 0 <= baseline - Fraction(255, 256) < Fraction(1, 4)
        ordinary = [v / (1 << j) for j, pair in enumerate(groups[level:]) for v in pair]
        for pos, h in enumerate(heights[:2*len(groups)]):
            group = pos // 2
            base = h / (1 << max(group-level, 0))
            assert 2 * base <= 1 and (base * (1 << height)).denominator == 1
        minimum = frozen + sum(ordinary) + unit * baseline
        reserve = unit / 4
        for sign in (-1, 0, 1):
            wanted = intercept + level * step + sign * step / MARGIN
            residual = wanted - minimum - reserve
            assert 0 <= residual <= sum(groups[level]) + target
            paid = minimum
            for value in ordinary:
                if value <= residual:
                    residual -= value
                    paid += value
            word = (residual + reserve) * (1 << digits) / unit
            assert word.denominator == 1 and 1 << (digits-2) <= word <= 3 << (digits-2)
            word = int(word)
            selected = 0
            for j in range(digits, 8, -1):
                if word % 2:
                    word -= units[j-1]
                    selected += units[j-1] << (digits-j)
                word //= 2
            assert 0 <= word < 256
            selected += word << (digits-8)
            assert paid + unit * Fraction(selected, 1 << digits) == wanted
        frozen += 2 * sum(groups[level])
    assert max_used == height
    print("87 levels; 261 exact mass reconstructions; maximum dyadic exponent 3360.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    main()

