"""Exact certificate for a nonrigid modulus with 196 squarefree odd factors."""
from fractions import Fraction
from hashlib import sha256
from math import isqrt, prod

from sympy import factorint, primerange

from ppi_mixed_slots import verify_templates
from ppi_partial_upper_pair import apply_partial_blocks
from ppi_rank_box_pair import odd_part
from ppi_weighted_fibers import bounds, tag_parameters


def prime_by_division(n):
    return n >= 2 and all(n % k for k in range(2, isqrt(n) + 1))


def add_squarefree_stars(used, count):
    """Choose fresh arms p,q,r satisfying (p-1)+(q-1)=r-1."""
    digest = sha256()
    triples = []
    candidates = list(primerange(2000, 4000))
    for i, p in enumerate(candidates):
        if p in used:
            continue
        for q in candidates[i+1:]:
            r = p + q - 1
            if q in used or r in used or not prime_by_division(r):
                continue
            assert p < q < r and prime_by_division(p) and prime_by_division(q)
            assert (p - 1) + (q - 1) == r - 1
            used.update((p, q, r))
            triples.append((p, q, r))
            digest.update(f"{p},{q},{r}\n".encode("ascii"))
            break
        if len(triples) == count:
            break
    assert len(triples) == count
    print(f"Squarefree damping triples: {count}; first {triples[0]}; last {triples[-1]}.")
    return digest.hexdigest()


def main():
    fixed = apply_partial_blocks(bounds(), [(1301, 1303), (179, 181)])
    axes = fixed[0]
    assert len(axes) == 52 and fixed[-2] == 48
    raw = [p for p in primerange(3, 1000) if p not in axes][:7]
    used = set(raw) | axes
    compensated = []
    for p in primerange(3, 30000):
        factors = factorint(odd_part(p-1))
        if p in used or set(factors) & used or any(e != 1 for e in factors.values()):
            continue
        compensated.append(p)
        used.add(p)
        used.update(factors)
        if len(compensated) == 18:
            break
    assert len(compensated) == 18
    data = tag_parameters(fixed, raw, compensated[4:], compensated[:4], 7, 3)
    tags, aux, tag_order, _, _, fine_weights = data
    aux_factors = factorint(aux)
    assert len(tags) == 25 and len(aux_factors) == 17
    assert all(e == 1 and prime_by_division(q) for q,e in aux_factors.items())
    assert set(aux_factors).isdisjoint(axes | set(tags))
    assert prod(aux_factors) == aux and prod(tags) == tag_order
    assert all(aux % odd_part(w) == 0 for w in fine_weights)
    total = verify_templates(fixed, *data, 700, 34, Fraction(3), head_count=320,
                             damping_builder=add_squarefree_stars)
    assert total == 196


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    main()

