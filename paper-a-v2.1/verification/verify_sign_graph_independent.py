#!/usr/bin/env python3
"""Independent exhaustive radial/ideal sign-graph audit; standard library only."""
from collections import defaultdict
from itertools import product
from math import gcd
from pathlib import Path
import json
import sys

def factor(n):
    ans=[]
    p=2
    while p*p<=n:
        if n%p==0:
            a=0
            while n%p==0:
                n//=p;a+=1
            ans.append((p,a))
        p+=1
    if n>1:ans.append((n,1))
    return ans

def divisors(n):
    return sorted(d for d in range(1,n+1) if n%d==0)

def mobius(n):
    fs=factor(n)
    if any(a>1 for _,a in fs):return 0
    return (-1)**len(fs)

def ramanujan(d,k):
    return sum(j*mobius(d//j) for j in divisors(gcd(d,k)))

def val(n,p):
    s=0
    while n%p==0:
        n//=p;s+=1
    return s

def components(n,ds,F):
    idx={c:i for i,c in enumerate(ds)}
    adj={c:[] for c,f in zip(ds,F) if f}
    for p,a in factor(n):
        for c in ds:
            s=val(c,p)
            if s<a and c in adj and c*p in adj:
                m=p**(a-s)
                if (2*F[idx[c]])%m:
                    adj[c].append(c*p);adj[c*p].append(c)
    comps=[];seen=set()
    for c in adj:
        if c in seen:continue
        todo=[c];seen.add(c);comp=[]
        while todo:
            x=todo.pop();comp.append(x)
            for y in adj[x]:
                if y not in seen:seen.add(y);todo.append(y)
        comps.append(comp)
    return comps

def autocorr(bits):
    n=len(bits)
    return tuple(sum(bits[x]*bits[(x+t)%n] for x in range(n)) for t in range(n))

def run(bound):
    totals=dict(moduli=bound,radial_profiles=0,subgroup_unions=0,direct_time_profiles=0,radial_count_checks=0,subgroup_count_checks=0)
    per=[]
    for n in range(1,bound+1):
        ds=divisors(n);L=len(ds)
        mat=[[ramanujan(d,n//c) for d in ds] for c in ds]
        downward=[sum(1<<j for j,e in enumerate(ds) if d%e==0) for d in ds]
        signatures=defaultdict(list);ideal_signatures=defaultdict(list)
        records=[]
        for mask in range(1<<L):
            selected=[j for j in range(L) if mask>>j&1]
            F=tuple(sum(row[j] for j in selected) for row in mat)
            sig=tuple(f*f for f in F)
            ideal=all(mask&downward[j]==downward[j] for j in selected)
            signatures[sig].append(mask)
            if ideal:ideal_signatures[sig].append(mask)
            records.append((mask,F,sig,ideal))
        ctimes={}
        for mask,F,sig,ideal in records:
            totals['radial_profiles']+=1
            if ideal:totals['subgroup_unions']+=1
            if mask==0:
                assert len(signatures[sig])==len(ideal_signatures[sig])==1
                continue
            comps=components(n,ds,F)
            ur=sum(1 not in comp for comp in comps)
            assert len(signatures[sig])==2**ur,(n,mask,'radial',ur,len(signatures[sig]))
            totals['radial_count_checks']+=1
            # More than a count: the actual component masks equal all signature mates.
            original_sign={c:(F[j]>0) for j,c in enumerate(ds) if F[j]}
            seen_masks=set()
            for mate in signatures[sig]:
                G=records[mate][1]
                flip=[]
                for comp in comps:
                    signs={original_sign[c] != (G[ds.index(c)]>0) for c in comp}
                    assert len(signs)==1
                    flip.append(signs.pop())
                assert not any(flip[k] for k,comp in enumerate(comps) if 1 in comp)
                seen_masks.add(tuple(flip))
            assert len(seen_masks)==2**ur
            if ideal:
                ud=sum(not any(c==1 or len(factor(c))==1 for c in comp) for comp in comps)
                assert len(ideal_signatures[sig])==2**ud,(n,mask,'ideal',ud,len(ideal_signatures[sig]))
                totals['subgroup_count_checks']+=1
            if n<=30:
                bits=[(mask>>ds.index(n//gcd(x,n)))&1 for x in range(n)]
                C=autocorr(bits)
                if sig in ctimes:assert ctimes[sig]==C
                else:ctimes[sig]=C
                totals['direct_time_profiles']+=1
        if n<=30:
            # Distinct Fourier signatures must also have different direct time data.
            assert len(set(ctimes.values()))==len(ctimes)
        per.append(dict(N=n,radial=1<<L,ideal=sum(r[3] for r in records),largest_radial_mate_family=max(map(len,signatures.values())),largest_ideal_mate_family=max(map(len,ideal_signatures.values()))))
    return dict(status='PASS',method='Exhaustive physical-order-shell enumeration with exact Ramanujan transforms; direct cyclic autocorrelation on every nonempty radial profile for N<=30.',totals=totals,per_modulus=per)

if __name__=='__main__':
    bound=int(sys.argv[1]) if len(sys.argv)>1 else 60
    result=run(bound)
    print(json.dumps(result,indent=2))
