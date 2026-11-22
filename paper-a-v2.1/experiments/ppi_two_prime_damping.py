"""Enumerate disjoint prime-power damping blocks by exact primality tests."""
from hashlib import sha256
from sympy import isprime, primerange

def main():
    used = set()
    digest = sha256()
    count = 0
    for p in primerange(2**21, 2**24):
        q = p * (p - 1) + 1
        if p in used or q in used or not isprime(q):
            continue
        assert q < 2**64
        used.update((p, q))
        digest.update(f"{p},{q}\n".encode("ascii"))
        count += 1
        if count == 12393:
            print(count, p, q, digest.hexdigest())
            return
    raise ValueError("insufficient blocks")

if __name__ == "__main__":
    main()
