"""Exact disjoint coarse/fine spectator certificate for the block-rank pair."""
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod
from sympy import factorint, primerange

from ppi_modular_pair import modular_bounds
from ppi_reallocated_pair import add_damping
from ppi_rank_box_pair import odd_part


def weight_counts(primes):
    result = [1]
    for p in primes:
        result = [x+y for x,y in zip(result+[0],[0]+[(p-1)*v for v in result])]
    return result


def main():
    bounds = modular_bounds()
    axes = bounds[0]
    smooth = {13,19,29}
    safe = [p for p in primerange(3,100000) if p not in axes|smooth
            and set(factorint(p-1)) <= smooth|{2}]
    raw = [p for p in primerange(3,1000) if p not in axes|smooth|set(safe)]
    total = verify_slots(bounds, safe[0], safe[1:18], raw[:16], 8, 4, 1100, 80, 128)
    assert total == 242


def verify_slots(bounds, guard, fine, gs, r, d, limit, t, margin):
    axes = bounds[0]
    tags = [guard]+fine+gs
    assert len(set(tags)) == len(tags) and axes.isdisjoint(tags)
    assert all(p > 2 and all(p % k for k in range(2,isqrt(p)+1)) for p in tags)
    m = len(gs)
    fine_weights = [prod(p-1 for p in t)*(guard-1) for t in combinations(fine,r)]
    aux = prod(odd_part(p-1) for p in [guard]+fine)
    assert gcd(aux,prod(axes)*guard*prod(fine)) == 1
    fine_factor = prod(p-1 for p in fine)
    wc = weight_counts(fine)
    K = guard*prod(fine)*prod(gs)
    raw_counts = weight_counts(gs)
    # Coarse slots miss the guard and exactly d ordinary axes.
    ordinary = sorted((fine_factor*prod(p-1 for p in t)
                       for t in combinations(gs,m-d)),reverse=True)
    n = (prod(fine)*sum(raw_counts[:m-d+1])-sum(ordinary)
         +(guard-1)*sum(wc[:r]))
    return verify_templates(bounds, tags, aux, K, n, ordinary, fine_weights, limit, t, margin)


def verify_templates(bounds, tags, aux, K, n, ordinary, fine_weights, limit, t, margin,
                     head_count=None, fine_units=None, prepared_groups=None, *,
                     damping_builder=add_damping):
    axes,E,dy,upper,lower,U,L,eu,el,fixed,_ = bounds
    axes = set(axes)
    # The four largest fixed primes are the lower modular helper.
    lower_helper = sorted(axes)[-4:]
    root = prod(q-1 for q in lower_helper)
    flat_exp = root.bit_length()+min(lower_helper).bit_length()+1
    alpha = Fraction((1 << flat_exp)//root,1 << flat_exp)
    eta = Fraction(1,min(lower_helper)-1)
    assert 0 <= 1-alpha*root <= eta and alpha*root <= 1
    base_bound = ((2*(1+U)+(1+U)*(1+L))/4
                  + (1+L)*eta)
    anchor_difference = 5+3*(1+upper)*(1+lower)
    M = min(v & -v for v in fine_weights)
    extra = max(((v & -v)//M).bit_length()-1 for v in fine_weights)
    assert len(ordinary) >= limit
    # Tail slots are first reached at positive binary depth. Their virtual
    # base heights may exceed one, provided every height actually used fits.
    head_count = limit if head_count is None else head_count
    assert 0 < head_count <= limit and head_count%2 == limit%2 == 0
    if prepared_groups is not None:
        target, ordered = prepared_groups
        assert len(ordered)*2 == limit and all(len(g) == 2 for g in ordered)
        assert all(target < 2*v <= 2*target for group in ordered for v in group)
        mean = sum(map(sum,ordered[:head_count//2]))/(head_count//2)
    else:
        target = 1 << (ordinary[head_count-1].bit_length()-2)
        ordered = []
        for start,end in ((0,head_count),(head_count,limit)):
            if start == end:
                continue
            masses = []
            for V in ordinary[start:end]:
                shift = (V-1).bit_length()-(target.bit_length()-1)
                v = Fraction(V,1 << shift) if shift >= 0 else Fraction(V << -shift)
                assert target < 2*v <= 2*target
                masses.append(v)
            masses.sort()
            groups = [(masses[i],masses[-1-i]) for i in range(len(masses)//2)]
            group_mean = sum(a+b for a,b in groups)/len(groups)
            if start == 0:
                mean = group_mean
            groups.sort(key=sum)
            drift = 0
            while groups:
                g = groups.pop(0 if drift>0 else -1)
                ordered.append(g)
                drift += sum(g)-group_mean
    s = max(0,(E*2**20//target).bit_length())
    step = Fraction(E*round(2*mean*2**s/E),2**s)
    delta = 2*step/K
    R = fixed+2*t
    if head_count < limit:
        first_depth = head_count//2-R
        assert first_depth >= 0 and 2*target <= ordinary[limit-1] << first_depth
    error = (base_bound+Fraction(3,2)*R*delta*(1+U)*(1+L))/2**t
    assert error + delta/(33*margin) < delta/margin, (t, float(error/delta))
    depth = 0
    while target > (M << depth) // 2:
        depth += 1
    assert R+depth+1 <= len(ordered)
    groups_used = ordered[:R+depth+1]
    # Find a common intercept leaving each residual in the greedy range.
    # Fine baseline at level l is (l+1)M-M/2^A; remove the tiny
    # term in the intercept when A has been determined.
    frozen = 0
    lows, highs = [], []
    for level in range(R+1):
        if fine_units is None:
            end = level+depth+1
            fine_estimate = (level+1)*M
        else:
            end = len(groups_used)
            fine_estimate = Fraction(5*M,4 << (R-level))
        base = (frozen + sum(sum(g)/2**j for j,g in
                 enumerate(groups_used[level:end])) + fine_estimate-level*step)
        lows.append(base+step/margin)
        highs.append(base+sum(groups_used[level])+target-step/margin)
        frozen += 2*sum(groups_used[level])
    lo, hi = max(lows), min(highs)
    assert lo < hi
    # Use the unadjusted midpoint. The omitted tiny term is
    # accounted for by the strict margin verified below.
    intercept = (lo+hi)/2
    c = 2*(n+intercept+R*step)/K
    assert 0 < c-delta*R <= c < 2 and 0 < delta < 1
    base_exp = (K*c).denominator.bit_length()-1
    assert (1 << base_exp) == (K*c).denominator
    blend_exp = ((33*margin)*E*anchor_difference//delta).bit_length()+1
    lam = Fraction(((1 << blend_exp)//E)*E,1 << blend_exp)
    assert (1-lam)*anchor_difference < delta/(33*margin) and 0 < lam <= 1
    A = max(blend_exp+t+max(base_exp+dy+flat_exp+3,s+dy,eu+el)+6,
            K.bit_length()+depth+5)+M.bit_length()
    if fine_units is None:
        assert A+R <= len(fine_weights) and Fraction(M,1 << A) < hi-intercept
    else:
        A = max(A, K.bit_length()+depth+R+5+M.bit_length())
        assert len(fine_units) == len(fine_weights) >= A
        assert all(w > 0 and w%2 for w in fine_units)
        units = fine_units[:A]
        prefix = next((j for j,w in enumerate(units) if w != 1), A)
        tail = sum(w << (A-j) for j,w in enumerate(units,1) if j > prefix)
        assert 4*tail < 1 << A
        baseline = Fraction((1 << A)-(1 << (A-prefix))+tail,1 << A)
    groups = groups_used
    # Verify exact endpoint decoding at every rank, including grid precision.
    frozen = 0
    for level in range(R+1):
        if fine_units is None:
            end = level+depth+1
            fine_min = (level+1)*M-Fraction(M,1 << A)
            reserve = 0
        else:
            end = len(groups)
            unit = Fraction(M,1 << (R-level))
            fine_min = unit*baseline
            reserve = unit/4
        digits = [v/2**j for j,g in enumerate(groups[level:end]) for v in g]
        base = frozen+sum(digits)+fine_min+reserve
        for sign in (-1,0,1):
            residual = intercept+level*step+sign*step/margin-base
            assert 0 <= residual <= sum(groups[level])+target
            for v in digits:
                if residual >= v:
                    residual -= v
            if fine_units is None:
                z = residual*(1 << A)/M
                assert z.denominator==1 and 0 <= z < 1 << A
            else:
                z = (residual+reserve)*(1 << A)/unit
                assert z.denominator == 1 and 1 << (A-2) <= z <= 3 << (A-2)
                numerator = int(z)
                for w in reversed(units[prefix:]):
                    if numerator%2:
                        numerator -= w
                    numerator //= 2
                assert 0 <= numerator < 1 << prefix
        frozen += 2*sum(groups[level])
    axes.update(tags+list(factorint(aux)))
    digest = damping_builder(axes,t)
    total = len(axes)
    print(f"axes={total}, t={t}, R={R}, limit={limit}, step-grid exponent={s}")
    print(f"A={A}, extra={extra}, depth={depth}, s={s}, damping digest={digest}")

    if fine_units is not None:
        print(f"Fine unit prefix={prefix}, available={len(fine_units)}, dyadic exponent={A+extra+R}")
    return total


if __name__ == "__main__":
    main()

