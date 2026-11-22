"""Check a six-axis two-mass certificate with exact totient counts."""
from math import gcd, prod


def check(primes, groups):
    points = [u for group in groups for u in group]
    assert len(groups) >= 2 and all(groups)
    assert len(set(points)) == len(points)
    assert len(set(map(sum, points))) == 1
    sizes = {sum(prod((p-1)*p**k for p, k in zip(primes, u))
                 for u in group) for group in groups}
    assert len(sizes) == 1
    return sizes.pop()


def main():
    first = [
        [(0,5,3),(0,6,2),(1,3,4)],
        [(1,0,7),(3,2,3),(6,0,2)],
        [(0,7,1),(0,8,0),(1,5,2),(1,6,1),(1,7,0),(2,4,2),(2,5,1),(2,6,0),(3,3,2)],
        [(1,2,5),(3,4,1),(6,2,0),(2,1,5),(4,3,1),(7,1,0),(3,0,5),(5,2,1),(8,0,0)],
    ]
    points = [(i,j,4-i-j) for i in range(5) for j in range(5-i)]
    second = [[u for j,u in enumerate(points) if mask >> j & 1]
              for mask in (0x6882, 0x1438)]
    a = check((3,5,7), first)
    b = check((17,23,47), second)
    assert (a,b) == (131418000,51166736192)
    assert gcd(a,b) == 16
    print(f"Verified masses {a} and {b}, gcd {gcd(a,b)}.")


if __name__ == "__main__":
    main()
