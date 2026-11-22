"""Verify coarse and fine shell allocations for the compressed PPI pair."""
from fractions import Fraction
from itertools import combinations, islice
from math import comb, gcd, prod
from sympy import factorint, primerange

from ppi_reallocated_pair import add_damping, donor_bounds
from ppi_reallocated_donors import UPPER, packet_bound
from ppi_reduced_lower import MASKS_1, SIGNED_1, PATH, pinned
from ppi_rank_box_pair import odd_part, blend

SHORT_PACKETS = [[[3615, 36119, 119183, 23135, 23447, 45903, 45979, 152639, 84543, 104751, 73047, 168795, 27483, 230813, 44446, 104381, 170735], 853], [[1477183, 1337919, 1172063, 322351, 942895, 840351, 770719, 63615, 193143, 647487, 3545583, 3546075, 3589099, 3304429, 3379966, 3053119, 3113469, 3838942, 1485310, 1633911, 745295], 4801]]


def check_sections(tags, r, R, A, extra, step, depth):
    """Check all rank templates in integer units of M/2^A."""
    vals = sorted(((p-1) & -(p-1)).bit_length()-1 for p in tags)
    b = sum(vals[:r])
    scale = 1 << (A-b)
    masses = (prod(p-1 for p in subset)*scale for subset in combinations(tags,r))
    step *= scale
    assert step.denominator == 1
    step = int(step)
    assert step % 128 == 0
    width = step//32

    def trim(mass, limit):
        exponent = max(0, mass.bit_length()-limit.bit_length())
        if mass > limit << exponent:
            exponent += 1
        assert exponent <= A+extra and mass % (1 << exponent) == 0
        weight = mass >> exponent
        assert weight <= limit < 2*weight
        return weight

    groups = [[trim(next(masses),width) for _ in range(4)]
              for _ in range(R+depth+1)]
    for mass in islice(masses,R+A):
        shell = mass//scale
        offset = ((shell & -shell).bit_length()-1)-b
        assert 0 <= offset <= extra
    frozen = permanent = 0
    for level in range(R+1):
        digits = []
        for j, group in enumerate(groups[level:level+depth+1]):
            assert all(v % (1 << j) == 0 for v in group)
            digits.extend(v >> j for v in group)
        # Fine weights are 2^(A-j) in these units.
        baseline = frozen+sum(digits)+(level+1)*(1 << A)-1
        center = (level+2)*step
        gap = center-baseline-permanent-width
        assert 0 <= gap <= 2*step
        for mass in islice(masses,10):
            if gap >= width//4:
                used = trim(mass,gap)
                permanent += used
                gap -= used
        assert 0 <= gap < width//4
        for offset in (-step//128,0,step//128):
            amount = center+offset-baseline-permanent
            initial = amount
            for digit in digits:
                if amount >= digit:
                    amount -= digit
            assert 0 <= amount < 1 << A
            assert baseline+permanent+(initial-amount)+amount == center+offset
        frozen += 2*sum(groups[level])


def main():
    axes, E, dyadic, upper, lower, U, L = donor_bounds()
    original = pinned(PATH)
    # Use shorter packets for the rational targets and moving corrections;
    # the original dyadic anchor and its bounds are unchanged.
    upper_bounds, ds = [], []
    for (hidden,fibers,_,scale,_),(masks,determinant) in zip(UPPER,SHORT_PACKETS):
        rows = [[i for i in range(len(fibers)) if m >> i & 1] for m in masks]
        d = original["donor"](hidden,fibers,rows,scale,determinant)
        upper_bounds.append(Fraction(packet_bound(original,fibers,rows,d),d))
        ds.append(d)
    U = max(upper_bounds)
    upper, eu = blend(upper,upper_bounds[0],ds[0],0)
    E *= odd_part(prod(ds))
    dyadic += max((d & -d).bit_length()-1 for d in ds)
    ld,lower_norm = pinned("scripts/ppi_lower_shield_pair.py")["inspect"](SIGNED_1,MASKS_1)
    lower, el = blend(lower,lower_norm,ld,25)
    select_sections(axes,E,dyadic,upper,lower,U,L,eu,el,79,3)


def select_sections(axes,E,dyadic,upper,lower,U,L,eu,el,fixed_rank,padding_rank):
    fixed_count = len(axes)
    tags = []
    excluded = prod(axes)
    for p in primerange(3, 100000):
        if gcd(p*(p-1),excluded) == 1:
            tags.append(p)
            excluded *= p
            if len(tags) == 34:
                break
    # Examine antichains within the same deterministic tag stream.
    candidates = []
    for size in range(19,len(tags)+1):
        ps = tags[:size]
        K = prod(ps)
        aux = prod(odd_part(p-1) for p in ps)
        aux_primes = set(factorint(aux))
        assert gcd(K*prod(axes),aux) == 1
        counts = [1]
        for p in ps:
            counts = [x+y for x,y in zip(counts+[0],[0]+[(p-1)*v for v in counts])]
        vals = sorted(((p-1)& -(p-1)).bit_length()-1 for p in ps)
        for r in range(size//2,size):
            capacity = comb(size,r)
            if capacity < 8000:
                continue
            minimum_mass = prod(p-1 for p in ps[:r])
            b = sum(vals[:r])
            M = 2**b
            extra = sum(vals[-r:])-b
            # step is measured in tag mass units, before division by K.
            s = (8*E//minimum_mass).bit_length()
            step = Fraction(E,2**s)
            assert minimum_mass >= 4*step
            if step <= 128*M:
                continue
            delta = 2*step/K
            n = sum(counts[:r])
            for t in range(1,1000):
                R = fixed_rank+padding_rank*t
                error = (4*(1+upper)*(1+lower)+6*R*delta*(1+U)*(1+L))/2**t
                if 128*error < delta:
                    break
            else:
                continue
            c = 2*(n+(R+2)*step)/K
            if not 0 < c < 2:
                continue
            # K times each free section density is dyadic. Coarse digits
            # have denominator at most 2^J and fine digits carry weight M.
            J = K.bit_length()+2
            A = max(t+s+eu+el+dyadic+4,J)+b+3
            W = step/32
            depth = 0
            while W > Fraction(M << depth,2):
                depth += 1
            demand = 4*(R+depth+1) + (R+A) + 10*(R+1)
            if demand <= capacity:
                candidates.append((fixed_count+2*t+size+len(aux_primes), t, size, r,
                                   A,extra,R,s,depth,demand,capacity))
    total,t,size,r,A,extra,R,s,depth,demand,capacity = min(candidates)
    tags = tags[:size]
    aux = prod(odd_part(p-1) for p in tags)
    check_sections(tags,r,R,A,extra,Fraction(E,2**s),depth)
    axes.update(tags)
    axes.update(factorint(aux))
    before = len(axes)
    digest = add_damping(axes,t)
    assert len(axes) == before+2*t
    assert total == len(axes)
    print(f"tags={size}, rank={r}, auxiliary={len(factorint(aux))}, t={t}, R={R}")
    print(f"s={s}, A={A}, extra={extra}, a={A+extra}, depth={depth}")
    print(f"slots={capacity}, demand={demand}, distinct odd primes={total}")
    print("damping SHA-256",digest)


if __name__ == "__main__":
    main()
