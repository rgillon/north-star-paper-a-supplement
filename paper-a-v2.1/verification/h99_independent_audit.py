#!/usr/bin/env python3
"""Independent finite checks for H99 and H100 rank truncation.
Uses only certificate DATA and already audited H100 donor exact rationals.
No imports from construction or verifier modules. Standard library only.
"""
import hashlib,json,sys
from fractions import Fraction as F
from itertools import product
from math import prod,isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parent
EXP=ROOT.parent/'experiments'
if len(sys.argv)>1:EXP=Path(sys.argv[1]).resolve()
DONOR=ROOT.parent/'evidence/h100_independent.json'
if len(sys.argv)>2:DONOR=Path(sys.argv[2]).resolve()
TAGS=(3,7,19,23,37,41,43,47,71,83,113,149)
POWERS=(5,3,2,2,2,2,1,1,1,1,1,1)
VP=(19,43,47,71,83,149)
VA=(4,3,3,3,3,2)
CP=(23,37,41)

def phi(p,a):return 1 if not a else (p-1)*p**(a-1)
def mass(ps,es):return prod(phi(p,a) for p,a in zip(ps,es))
def two(n):return (n&-n).bit_length()-1
def dbits(x):
    x=F(x);d=x.denominator
    assert d&(d-1)==0
    return d.bit_length()-1

def cp(n):
    z=1
    for p in CP:
        while n%p==0:n//=p;z*=p
    return z

def norm(w,H):
    b=F(1)
    while w*b>H:b/=2
    while 2*w*b<=H:b*=2
    assert H/2<w*b<=H
    return b

def truncation(tags,powers,guards,selected,R,A,H,s,b=None):
    h=83349;B=prod(p**a for p,a in zip(VP,VA));Q=prod(p**a for p,a in zip(tags,powers))
    K=Q//h*B
    orders=list(product(*(range(a+1) for a in powers)))
    weight=[mass(tags,v) for v in orders]
    assert sum(weight)==Q
    pos={v:j for j,v in enumerate(selected)}
    assert len(pos)==len(selected)
    g=tuple(a if p in guards else 0 for p,a in zip(tags,powers))
    bg={i for i,v in enumerate(orders) if i not in pos and any(v[j]<powers[j] for j,p in enumerate(tags) if p in guards)}
    index={v:i for i,v in enumerate(orders)}
    for i,v in enumerate(orders):
        if v!=g and all(x<=y for x,y in zip(v,g)):assert i in bg
        if i in pos:assert any(v[j]<powers[j] for j,p in enumerate(tags) if p in guards)
        if i not in bg and i not in pos:continue
        for j,a in enumerate(v):
            if not a:continue
            u=v[:j]+(a-1,)+v[j+1:];k=index[u]
            if i in bg:assert k in bg
            elif k not in bg:
                assert k in pos and pos[k]<pos[i]
                assert weight[i]>=2*weight[k]
    assert all(weight[i]%h==0 for i in selected)
    bgmass=sum(weight[i] for i in bg);assert bgmass%h==0
    n=bgmass//h*B
    weights=[weight[i]//h*B for i in selected]
    H0=min(weights[:174 if len(tags)==13 else 214])*99//100
    betas=[norm(w,H0) for w in weights]
    zs=[w*z for w,z in zip(weights,betas)]
    pairs=list(zip(zs[::2],zs[1::2]))[:R+160]
    assert len(pairs)==R+160
    gm=mass(tags,g)
    fine=[]
    for v in product(*(range(a+1) for a in VA)):
        V=gm*mass(VP,v);d=cp(V)
        assert 34891%d==0
        u=(V>>two(V))//(h*d)
        assert u*h*d==(V>>two(V)) and u>0 and u%2
        fine.append((u,sum(v),v,V,d))
    fine.sort()
    fp={a[2]:i for i,a in enumerate(fine)}
    for i,(_,_,v,V,d) in enumerate(fine):
        for k,a in enumerate(v):
            if a:
                j=fp[v[:k]+(a-1,)+v[k+1:]]
                assert j<i and d%fine[j][4]==0 and V%fine[j][3]==0
    M=1<<min(two(a[3]) for a in fine)
    extra=[two(a[3])-two(M) for a in fine]
    assert min(extra)==0 and max(extra)==7
    units=[a[0] for a in fine];assert units[:8]==[1]*8
    assert H0/F(2**159)<=F(M,2)
    lows=[];highs=[];frozen=F(0)
    margin=F(27,10) if len(tags)==13 else F(7,2)
    for l in range(R+1):
        minimum=frozen+sum(sum(g)/2**j for j,g in enumerate(pairs[l:]))
        base=minimum+F(5*M,4*2**(R-l))-l*s
        lows.append(base+s/margin)
        highs.append(base+sum(pairs[l])+H0-s/margin)
        frozen+=2*sum(pairs[l])
    assert max(lows)<min(highs)
    mid=(max(lows)+min(highs))/2
    if b is not None:assert b==mid
    b=mid
    frozen=F(0);maxdepth=0;replays=0;covercount=0
    for l in range(R+1):
        off=R-l;D=A-off;unit=F(M,2**off)
        assert D>=8
        baseline=sum(F(u,2**j) for j,u in enumerate(units[:D],1))
        assert 0<=baseline-F(255,256)<F(1,4)
        for j in range(1,D+1):maxdepth=max(maxdepth,off+j+extra[j-1])
        for i,beta in enumerate(betas[:2*len(pairs)]):
            k=i//2
            ht=2*beta if k<l else beta/F(2**(k-l))
            assert ht<=1 and (k<l or 2*ht<=1)
            assert dbits(ht)<=H
        digits=[v/F(2**j) for j,g in enumerate(pairs[l:]) for v in g]
        minmass=frozen+sum(digits)+unit*baseline
        endpoints=[b+l*s-s/margin,b+l*s+s/margin]
        for wanted in endpoints:
            residue=wanted-minmass-unit/4
            assert 0<=residue<=sum(pairs[l])+H0
            covercount+=1
        # Select first and last valid absolute-grid targets, plus central grid targets.
        # The interval endpoints need not themselves be dyadic for H99.
        tick=F(M,2**A)
        lo=(endpoints[0]/tick).__ceil__();hi=(endpoints[1]/tick).__floor__()
        assert lo<=hi
        for integer in sorted(set((lo,lo+1,(lo+hi)//2,hi-1,hi))):
            if not lo<=integer<=hi:continue
            wanted=integer*tick
            residue=wanted-minmass-unit/4;paid=minmass
            for amount in digits:
                if amount<=residue:residue-=amount;paid+=amount
            word=(residue+unit/4)*2**D/unit
            assert word.denominator==1 and 2**(D-2)<=word<=3*2**(D-2)
            word=int(word);selectedmass=F(0)
            for j in range(D,8,-1):
                if word%2:
                    word-=units[j-1];selectedmass+=F(units[j-1],2**j)
                word//=2
            assert 0<=word<256
            selectedmass+=F(word,256)
            assert paid+unit*selectedmass==wanted
            replays+=1
        frozen+=2*sum(pairs[l])
    assert maxdepth==H
    c=2*(n+b+R*s)/K;delta=2*s/K
    assert 0<c-R*delta<c<2
    return dict(K=K,n=n,H0=H0,s=s,b=b,c=c,delta=delta,M=M,levels=R+1,precision=A,height=H,background=len(bg),fine=len(fine),endpoint_checks=covercount,independent_replays=replays,interval_width=(min(highs)-max(lows))/s)

def serial(x):
    if isinstance(x,F):return {'numerator':str(x.numerator),'denominator':str(x.denominator),'decimal':float(x)}
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    return x

manifest=json.loads((EXP/'ppi_100_certificate.json').read_text())
donor=json.loads(DONOR.read_text());assert donor['result']=='PASS'
E=int(donor['donor']['E'])
base_s=F(int(manifest['step_numerator']),2**manifest['step_binary_exponent'])
base_b=F(int(manifest['intercept_numerator']),2**manifest['intercept_binary_exponent'])
x100=truncation(TAGS,POWERS,(3,7,113),json.loads((EXP/'ppi_100_slot_order.json').read_text()),86,3353,3360,base_s,base_b)
tags=TAGS+(65537,);powers=POWERS+(1,)
assert all(65537%p for p in range(2,isqrt(65537)+1))
K=prod(p**a for p,a in zip(tags,powers))//83349*prod(p**a for p,a in zip(VP,VA))
s=E*F(round(F(130352847*K*2**192,10**12*E)),2**192)
x99=truncation(tags,powers,(3,7,113,65537),json.loads((EXP/'ppi_99_slot_order.json').read_text()),84,3401,3408,s)
assert x99['K']==K
U=F(**{k:int(v) for k,v in donor['exact_rationals']['upper_norm'].items()})
L=F(**{k:int(v) for k,v in donor['exact_rationals']['lower_norm'].items()})
primepowers={int(p):a for p,a in manifest['odd_prime_powers'].items()}
last=2101327;q=last*(last-1)+1
assert primepowers.pop(last)==2 and primepowers.pop(q)==1 and 65537 not in primepowers
primepowers[65537]=1;assert len(primepowers)==99
lowerhelpers=sorted(primepowers)[-4:]
flat=prod(p-1 for p in lowerhelpers).bit_length()+min(lowerhelpers).bit_length()+1
eta=F(1,min(lowerhelpers)-1)
alpha=F(2**flat//prod(p-1 for p in lowerhelpers),2**flat)
assert 0<=1-alpha*prod(p-1 for p in lowerhelpers)<=eta
blend=F(E*(2**400//E),2**400)
error=((1+U)*(3+L)/4+(1+L)*eta+F(3,2)*84*x99['delta']*(1+U)*(1+L))/2**19
# A rigorous upper bound slightly larger than the exact source anchor norm.
anchor=5+3*(1+602953959433)*(1+(L+1))
error+=(1-blend)*anchor
assert error/x99['delta']<F(339,1000)<F(10,27)
dy=donor['donor']['dy'];eu=donor['donor']['eu'];el=donor['donor']['el']
profile=400+19+max(dbits(K*x99['c'])+dy+flat+3,max(dbits(s),dbits(K*x99['delta']))+dy,eu+el)+6+x99['M'].bit_length()
ordinary=K.bit_length()+159+84+5+x99['M'].bit_length()
assert max(profile,ordinary)==3401
x99.update(error_over_delta=error/x99['delta'],allowed_error=F(10,27),precision_profile=profile,precision_ordinary=ordinary,flat_helper_bits=flat,odd_distinct=len(primepowers),bit_length=(prod(p**a for p,a in primepowers.items())<<3408).bit_length())
out={'result':'PASS','scope':'Independent changed arithmetic for H100 rank truncation and H99; previously audited H100 donor certificate imported as exact data, no construction/verifier modules imported.','H100_rank':x100,'H99':x99,'inputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [DONOR,EXP/'ppi_100_certificate.json',EXP/'ppi_100_slot_order.json',EXP/'ppi_99_slot_order.json']}}
print(json.dumps(serial(out),indent=2))
