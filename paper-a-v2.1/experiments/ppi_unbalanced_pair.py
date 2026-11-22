"""Exact Level 1 certificate using two-prime partial upper interpolants."""
from fractions import Fraction
from math import isqrt

from ppi_partial_upper_pair import apply_partial_blocks
from ppi_weighted_digits import spectator_data, verify_templates
from ppi_weighted_fibers import donor_bounds


HIDDEN = [5, 11, 13, 17]
HIDDEN_POWERS = [2, 2, 1, 1]
FIBERS = [[-89, 2, 1, 34359738368],
 [-103, 2, 1, 38654853120],
 [67, 4, 2, 67721379840],
 [241, 1, 1, 67645734912],
 [-73, 3, 2, 38730334208],
 [-257, 1, 2, 68585517056],
 [53, 5, 1, 67777855488],
 [-137, 2, 2, 68717527040],
 [271, 1, 1, 68660903936],
 [-131, 2, 1, 67777855488],
 [-97, 2, 1, 38654705664],
 [127, 2, 1, 67645882368],
 [-251, 1, 2, 67645865984],
 [-263, 1, 1, 68717379584],
 [139, 2, 1, 68712349696],
 [193, 1, 1, 38654705664]]
MASKS = [2656,
 960,
 16936,
 6528,
 11008,
 6852,
 32986,
 41062,
 43141,
 37652,
 57396,
 2980,
 52384,
 51746,
 19539]
PAIRS = [[776561, 776563], [881, 883], [29, 31]]
RAW = [3, 7, 23, 37, 41, 71]
GUARDS = [19, 43, 47, 83, 109, 113]
FINE = [149, 197, 337, 379, 433, 449, 487]


def bounds():
    base = donor_bounds(HIDDEN,FIBERS,MASKS,2400,193,249391 << 32,
                        (304,328,170,462),HIDDEN_POWERS)
    assert base[5] < 602953959432
    result = apply_partial_blocks(base,PAIRS)
    assert 30 < result[5] < 31
    assert len(result[0]) == 48 and result[-2] == 46
    return result


def main():
    args = bounds()
    assert all(p > 2 and all(p % k for k in range(2,isqrt(p)+1))
               for p in RAW+GUARDS+FINE)
    data = spectator_data(args[0],RAW,GUARDS,FINE,[2,2,2,2,2,1],3,12,3)
    total = verify_templates(args,*data[:-1],408,23,Fraction(5,2),
                             head_count=190,fine_units=data[-1])
    assert total == 113


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
