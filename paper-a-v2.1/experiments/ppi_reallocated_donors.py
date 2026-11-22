"""Reproduce the 79-axis fixed donor check using exact integer arithmetic."""
from math import gcd, prod
from ppi_reduced_lower import SIGNED_1, SIGNED_2, MASKS_1, MASKS_2, pinned, PATH


# Hidden primes, signed quotient primes with upset masks, conductor masks,
# mass scale, and signed incidence minor.
UPPER = [
    ((3, 7, 13, 19),
     [[433, 32768], [457, 32896], [541, 43008], [547, 43520], [571, 43648], [761, 57568], [739, 52416], [709, 51328], [577, 43688], [-827, 64160], [-859, 65484], [-863, 65516], [-787, 60032], [-823, 64128], [-811, 64512], [-857, 65248], [-809, 61088], [-853, 65216]],
     [52987, 248190, 84543, 85885, 27483, 250333, 81374, 125374, 53223, 152639, 58199, 63966, 167262, 175083, 109815, 154590, 23447],
     6, -853),
    ((5, 11, 17, 29, 31),
     [[3361, 2147483648], [3571, 2281701376], [3593, 2155905024], [3691, 2290089984], [3803, 2290122752], [3907, 2818572288], [4019, 2818605056], [4201, 3221225472], [4027, 2826960896], [3697, 2684354560], [5051, 3970465792], [-5147, 4206936256], [-5179, 4270376960], [-5189, 4277698560], [-5197, 4289382400], [-5233, 4293577864], [-5107, 4204824704], [-5153, 4240998400], [-5081, 3974660096], [-5227, 4294623360], [-5009, 3904929792], [-4801, 3766534272]],
     [684159, 3589099, 3594206, 2844223, 3751925, 3379966, 1477183, 2223679, 3838942, 745295, 3304429, 3113469, 647487, 3546075, 2225134, 193143, 1931260, 1485310, 1337919, 745271, 718015],
     160, 9602),
]


def packet_bound(upper, fibers, rows, denominator):
    primes = [abs(p) for p, _ in fibers]
    matrix = [[int(i in row) for row in rows] + [-denominator // p]
              for i, p in enumerate(primes)]
    _, triangular = upper["echelon"](matrix[:-1])
    coeff = [0] * len(rows)
    for i in reversed(range(len(rows))):
        rhs = triangular[i][-1] - sum(triangular[i][j]*coeff[j]
                                      for j in range(i+1, len(rows)))
        coeff[i], rem = divmod(rhs, triangular[i][i])
        assert rem == 0
    assert all(sum(x*y for x,y in zip(row[:-1],coeff)) == row[-1] for row in matrix)
    return abs(sum(coeff)) + sum(abs(a)*prod(primes[i]-1 for i in row)
                                 for a,row in zip(coeff,rows))


def main():
    upper = pinned(PATH)
    axes = list(map(abs, SIGNED_1 + SIGNED_2))
    denominators, bounds = [], []
    for hidden, fibers, masks, scale, det in UPPER:
        rows = [[i for i in range(len(fibers)) if mask >> i & 1] for mask in masks]
        assert all(0 < mask < 1 << len(fibers) for mask in masks)
        d = upper["donor"](hidden, fibers, rows, scale, det)
        denominators.append(d)
        bounds.append(packet_bound(upper, fibers, rows, d))
        axes.extend(hidden)
        axes.extend(abs(p) for p, _ in fibers)
    assert len(axes) == len(set(axes)) == 79
    d1, d2 = denominators
    assert gcd(d1,d2) == 1
    k = d1 * pow(d1,-1,d2)
    if 2*k > d1*d2:
        k -= d1*d2
    assert k % d1 == 0 and (1-k) % d2 == 0
    bound = abs(k//d1)*bounds[0] + abs((1-k)//d2)*bounds[1]
    lower = pinned("scripts/ppi_lower_shield_pair.py")
    (d1,b1),(d2,b2) = [lower["inspect"](c,m) for c,m in
                         [(SIGNED_1,MASKS_1),(SIGNED_2,MASKS_2)]]
    u,v = d1//(d1&-d1), d2//(d2&-d2)
    assert gcd(u,v) == 1 and max(d1&-d1,d2&-d2) == 2**25
    k = u*pow(u,-1,v)
    assert abs(k)*b1+abs(1-k)*b2 < 2**422
    print("79 distinct axes; upper sizes 22 and 27; lower sizes 13 and 17.")
    print("Integer upper witness sup norm < 2^" + str(bound.bit_length()))
    print("Lower denominator divides 2^25, lower sup norm < 2^422.")


if __name__ == "__main__":
    main()
