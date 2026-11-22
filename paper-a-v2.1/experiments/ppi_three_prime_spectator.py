"""Small exact demonstration of equal-mass antichains on three prime axes."""
from itertools import combinations, product
from math import prod

PRIMES = (3, 5, 7)


def layer(existing):
    for degree in range(1, 17):
        points = [(i, j, degree-i-j) for i in range(degree+1)
                  for j in range(degree-i+1)]
        weights = [prod(p**k for p, k in zip(PRIMES, u)) for u in points]
        seen = {}
        for size in range(1, 4):
            for subset in combinations(range(len(points)), size):
                weight = sum(weights[i] for i in subset)
                if weight in seen:
                    left = set(seen[weight]) - set(subset)
                    right = set(subset) - set(seen[weight])
                    digits = [points[i] for i in left | right]
                    sums = {tuple(x+y for x,y in zip(u,v))
                            for u,v in product(existing, digits)}
                    if len(sums) == len(existing)*len(digits):
                        return [[points[i] for i in sorted(side)] for side in (left, right)]
                seen[weight] = subset
    raise ValueError("bounded demonstration found no relation")


def main():
    collections = [[(0, 0, 0)]]
    for _ in range(2):
        sides = layer([u for group in collections for u in group])
        collections = [[tuple(x+y for x, y in zip(u, v))
                        for u, v in product(group, side)]
                       for group, side in product(collections, sides)]
    points = [u for group in collections for u in group]
    assert len(points) == len(set(points))
    assert len(set(map(sum, points))) == 1
    masses = [sum(prod((p-1)*p**k for p, k in zip(PRIMES, u))
                  for u in group) for group in collections]
    assert len(set(masses)) == 1
    print(len(collections), masses[0], max(max(u) for u in points)+1)
    print(collections)


if __name__ == "__main__":
    main()
