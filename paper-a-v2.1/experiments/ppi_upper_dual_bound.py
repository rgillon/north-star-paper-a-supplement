"""Full zero spectrum and a dual norm certificate for the PPI-WF1 upper block."""
from fractions import Fraction
from math import lcm, prod

from ppi_weighted_fibers import FIBERS, HIDDEN, LOWER, bounds


DUAL = (738,604,3020,4161,1007,7450,-268,5772,-8188,-5503,0,9597,-201,
        7181,-10000,-7987,5973,-4631,-470,-1946)


def subset_sums(values):
    result = [0]
    for value in values:
        result += [s+value for s in result]
    return result


def zero_masks():
    primes = [abs(p) for p,_,_ in FIBERS]
    common = lcm(*(p-1 for p in primes))
    result = []
    split = len(primes)//2
    for v in range(1 << len(HIDDEN)):
        hidden = [prod((-1 if v >> j & 1 else q-1)
                       for j,q in enumerate(HIDDEN) if shell >> j & 1)
                  for shell in range(1 << len(HIDDEN))]
        coefficients = []
        for signed,_,mask in FIBERS:
            p = abs(signed)
            raw = sum(w for j,w in enumerate(hidden) if mask >> j & 1)
            coefficients.append((signed//p)*raw*(common//(p-1)))
        target = sum(coefficients)
        weights = [p*a for p,a in zip(primes,coefficients)]
        left = {}
        for t,w in enumerate(subset_sums(weights[:split])):
            left.setdefault(w,[]).append(t)
        for t,w in enumerate(subset_sums(weights[split:])):
            result.extend((v,u | (t << split)) for u in left.get(target-w,[]))
    return result


def upper_lower_bound():
    zero = zero_masks()
    assert [(v,t) for v,t in zero if v] == [(16,912639)]
    quotient = [t for v,t in zero if not v]
    assert len(quotient) == len(set(quotient)) == 1469
    top = [h for h in range(64) if any(m >> h & 1 for _,_,m in FIBERS)]
    assert {bool(h & 16) for h in top} == {True,False}
    primes = [abs(p) for p,_,_ in FIBERS]
    # Transpose of packet coefficient extraction, with denominator prod(primes).
    values = [0]*(1 << len(primes))
    for t in quotient:
        values[t] = sum(y for i,y in enumerate(DUAL) if t >> i & 1)
    for i,p in enumerate(primes):
        stride = 1 << i
        for start in range(0,len(values),2*stride):
            for j in range(start,start+stride):
                root,detail = values[j],values[j+stride]
                values[j] = (p-1)*root+detail
                values[j+stride] = root-detail
    denominator = sum(abs(x) for x in values)
    numerator = abs(sum(y*(prod(primes)//p) for y,p in zip(DUAL,primes)))
    result = Fraction(numerator,denominator)
    assert 23518 < result < 23519
    return result


def main():
    lower_bound = upper_lower_bound()
    axes = bounds()[0]
    lower_helper = sorted(axes)[-4:]
    upper_helper = sorted(axes)[-8:-4]
    l = list(map(abs,LOWER[0]))
    lower_cost = Fraction(2*sum(p-1 for p in l),max(l)-min(l))
    upper_cost = Fraction(2*max(upper_helper),max(upper_helper)-min(upper_helper))
    z = prod(q-1 for q in lower_helper)
    vmin = -(z//(min(lower_helper)-1))
    vmax = -(z//(max(lower_helper)-1))
    helper_cost = Fraction(2*(1-vmin),vmax-vmin)
    cost = lower_cost*upper_cost*helper_cost
    assert cost < 112
    for t in range(14):
        assert (lower_bound-2)/(cost*2**t) > Fraction(5,16)*Fraction(2,44+2*t)
    print("1469 quotient zero masks, one mixed zero mask. Upper witness norm > 23518.")
    print("Contraction cost < 112 before damping. All 14 rank-window obstructions pass.")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run this certificate with Python assertions enabled.")
    main()
