"""Check the joint prime budget and both lower interpolants with exact arithmetic."""
from math import gcd
from ppi_reduced_upper import SOURCE, PATH
import subprocess

UPPER_INDICES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 18, 19, 20, 89, 99, 96, 98, 13, 14, 21, 22, 24, 29, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 92, 97]
UPPER_MASKS = [1474115584, 2216571999232, 137673068800, 207319991296, 1103999275024, 34474626048, 13194296606720, 171928584193, 276083706880, 8798374592648, 275010708480, 86052208640, 2207847942160, 6734600750080, 86066856576, 13744071286784, 70187247616, 1138468175872, 2199538897408, 138227302528, 1125365072000, 8798364139776, 429170688, 9895839404048, 10317761536, 2199666836480, 30685495296, 141863321856, 283662786560, 4415417139200, 19557779968, 4947936175104, 1172626510848, 276077093888, 1237075788032, 13211354923008, 1104066515456, 137671610384, 17591917871103, 139686838788, 4426828064, 1103923773512, 2199206625282]

SIGNED_1 = [23, 37, 43, 73, 79, 83, 109, 167, -41, -59, -157, -173, -181]
MASKS_1 = [1579, 2590, 2977, 3238, 3757, 3763, 4397, 4744, 5531, 5575, 6128, 7416]
SIGNED_2 = [179, 359, 839, 389, 479, 257, 797, 647, 509, 557, 1193, 503, 227, -3299, -1367, -977, -1283]
MASKS_2 = [8319, 8637, 10979, 12623, 14584, 14818, 27838, 28013, 28115, 30704, 32534, 41470, 44403, 75582, 79191, 121771]

def pinned(path):
    namespace = {"__name__": "pinned_verifier"}
    exec(compile(subprocess.check_output(["git", "show", f"{SOURCE}:{path}"]),
                 path, "exec"), namespace)
    return namespace


def main():
    upper = pinned(PATH)
    lower = pinned("scripts/ppi_lower_shield_pair.py")
    fibers = [(abs(upper["SECOND_FIBERS"][i][0]) * (1 if j < 22 else -1),
               upper["SECOND_FIBERS"][i][1]) for j,i in enumerate(UPPER_INDICES)]
    rows = [[j for j in range(44) if m >> j & 1] for m in UPPER_MASKS]
    first = upper["donor"](upper["FIRST_HIDDEN"], upper["FIRST_FIBERS"],
                           upper["FIRST_ZERO_ROWS"], 128, 151)
    second = upper["donor"](upper["SECOND_HIDDEN"], fibers, rows, 12, -869371)
    assert gcd(first, second) == 1
    axes = list(upper["FIRST_HIDDEN"] + upper["SECOND_HIDDEN"])
    axes.extend(abs(p) for p, _ in upper["FIRST_FIBERS"])
    axes.extend(abs(upper["SECOND_FIBERS"][i][0]) for i in UPPER_INDICES)
    axes.extend(abs(p) for p in SIGNED_1 + SIGNED_2)
    assert len(axes) == len(set(axes)) == 98
    (d1,b1), (d2,b2) = (lower["inspect"](c,m) for c,m in
                           [(SIGNED_1,MASKS_1),(SIGNED_2,MASKS_2)])
    u, v = d1 // (d1 & -d1), d2 // (d2 & -d2)
    assert gcd(u,v) == 1
    k = u * pow(u,-1,v)
    assert k % u == 0 and (1-k) % v == 0
    assert max(d1 & -d1,d2 & -d2) == 2**25
    bound = abs(k)*b1 + abs(1-k)*b2
    assert bound < 2**422
    print("98 disjoint odd primes; combined lower denominator divides 2^25.")
    print("Combined lower sup norm is less than 2^422.")


if __name__ == "__main__":
    main()
