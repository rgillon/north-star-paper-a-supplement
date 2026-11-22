"""Independent reconstruction of the last two-squared-axis population obstruction.

This check does not import the exhaustive grid verifier. It validates the
published two masks, their downset-difference realization, population
polynomial and largest-prime residue. It does not replace the exhaustive
enumeration of the remaining patterns.
"""
from itertools import product
from pathlib import Path
import json
import sympy as sp

p,q,r,s = sp.symbols('p q r s')
axes=(p,q,r,s)
cells=list(product(range(3),range(3),range(2),range(2)))
positive=18258070720
negative=9296846848
P={v for i,v in enumerate(cells) if positive >> i & 1}
Q={v for i,v in enumerate(cells) if negative >> i & 1}
assert P and Q and not P & Q

def below(u,v):
    return all(x<=y for x,y in zip(u,v))

D={v for v in cells if any(below(v,w) for w in P|Q)}
I={v for v in D if not any(below(w,v) for w in Q)}
J={v for v in D if not any(below(w,v) for w in P)}
assert I-J == P and J-I == Q
assert all(not below(u,v) or u in I for u in cells for v in I)
assert all(not below(u,v) or u in J for u in cells for v in J)

def population(v):
    return sp.prod(1 if k==0 else t**(k-1)*(t-1)
                   for t,k in zip(axes,v))

poly=sp.expand(sum(population(v) for v in P)-sum(population(v) for v in Q))
claimed=p**2*q**2*(r-s)-p**2*r*s+p**2*s+q**2*r*s-q**2*r
assert sp.expand(poly-claimed)==0
residue=sp.factor(poly.subs(s,0))
assert sp.expand(residue-q**2*r*(p-1)*(p+1))==0
report={'positive_mask':positive,'negative_mask':negative,
        'positive_shells':sorted(P),'negative_shells':sorted(Q),
        'positive_ideal_cells':len(I),'negative_ideal_cells':len(J),
        'population_polynomial':str(poly),'residue_mod_largest_prime':str(residue),
        'passed':True,'scope':'One last-pattern symbolic check, not exhaustive enumeration'}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
