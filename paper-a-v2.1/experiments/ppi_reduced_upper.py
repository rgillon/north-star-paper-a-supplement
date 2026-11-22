"""Run with Python 3 in the checkout containing the pinned PPI source commit."""
import math
import subprocess

SOURCE = "c8dbda833f437371efc12b9c1205006bef6cebac"
PATH = "scripts/ppi_upper_coprime_donors.py"
INDICES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 18, 19, 20, 91, 89, 99, 96, 98, 13, 14, 21, 22, 24, 28, 29, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 90, 92, 97]
CONDUCTORS = [35188728301568, 17592786617360, 13194969482240, 10995462086656, 8796728066050, 13194322771969, 17601132011520, 9886778368, 138451615808, 5057286656, 8796527396896, 43980815646720, 8846575620, 2375617536, 2886393398272, 17729971863552, 137668660224, 35287792910336, 5497680592912, 17867688427520, 2474523795456, 567114744832, 8801080999936, 343810515968, 69422286864, 38483216859136, 4402967257088, 18692013072384, 983318560, 69098717184, 5532528967680, 17592548427776, 829104325632, 4956622880, 36284501893120, 19793964728320, 17609585885440, 44530251106304, 35185021683840, 21990444138752, 35218852184096, 35734768517120, 70367670697983, 275008716932, 2346975304]

def main():
    source = subprocess.check_output(["git", "show", f"{SOURCE}:{PATH}"])
    donor = {"__name__": "pinned_donor"}
    exec(compile(source, PATH, "exec"), donor)
    pool = donor["SECOND_FIBERS"]
    fibers = []
    for j, i in enumerate(INDICES):
        p, bits = pool[i]
        fibers.append((abs(p) if j < 23 else -abs(p), bits))
    rows = [[j for j in range(46) if mask >> j & 1] for mask in CONDUCTORS]
    first = donor["donor"](donor["FIRST_HIDDEN"], donor["FIRST_FIBERS"],
                           donor["FIRST_ZERO_ROWS"], 128, 151)
    second = donor["donor"](donor["SECOND_HIDDEN"], fibers, rows, 12, -4346855)
    assert math.gcd(first, second) == 1
    primes = sorted(abs(p) for p, _ in pool)
    gaps = [sum(primes[:m]) - sum(primes[-(m-2):]) for m in range(3, 18)]
    assert min(gaps) == 23628
    print("Second donor: 46 quotient + 5 hidden. Fixed donors: 112.")
    print("Minor -4346855; denominator gcd 1; minimum extremal gap 23628.")


if __name__ == "__main__":
    main()
