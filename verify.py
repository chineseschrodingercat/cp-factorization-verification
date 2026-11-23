"""Exact four-bit checks with optional SymPy verification of the arithmetic examples."""
import argparse
from fractions import Fraction as Q
from itertools import combinations, product
from math import isqrt
from collections import defaultdict

WORDS = tuple(product((0, 1), repeat=4))
def parity(s):
    return (-1)**sum(s)
def readout(s):
    return s[0]*s[1] + s[2]*s[3]
def dot(a,b):
    return sum(x*y for x,y in zip(a,b))
def outer(u,v):
    return tuple(tuple(a*b for b in v) for a in u)
def sum_matrices(mats):
    return tuple(tuple(sum(row) for row in zip(*rs)) for rs in zip(*mats))
def rank(rows):
    rows = [list(map(Q,row)) for row in rows]
    r=0
    for j in range(len(rows[0])):
        pivot=next((i for i in range(r,len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        for i in range(r+1,len(rows)):
            if rows[i][j]:
                t=rows[i][j]/rows[r][j]
                rows[i]=[a-t*b for a,b in zip(rows[i],rows[r])]
        r+=1
        if r==len(rows):
            break
    return r
def local_test():
    bits=[[s[i] for s in WORDS] for i in range(4)]
    moment=[[1]*16,*bits]
    moment += [[x*y for x,y in zip(a,b)] for a,b in combinations(bits,2)]
    f=list(map(readout,WORDS)); chi=list(map(parity,WORDS))
    moment += [[x*y for x,y in zip(a,f)] for a in bits]
    assert rank(moment)==15 and all(dot(c,chi)==0 for c in moment)
    square=[x*x for x in f]
    assert rank([*moment,square])==16 and sum(square)==10 and dot(square,chi)==2
    print("Local moment ranks 15 -> 16; parity kernel and slope 2 verified.")

def rank_one_squares(weight):
    """Nonnegative rational squares summing to weight."""
    weight=Q(weight)
    num=weight.numerator*weight.denominator
    j=0
    while num:
        if num & 1:
            yield Q(1<<(j//2), weight.denominator)
            if j%2:
                yield Q(1<<(j//2), weight.denominator)
        num>>=1
        j+=1

def model(B):
    p,q=len(B),len(B[0]); m=p+q
    sides=(tuple(range(p)),tuple(range(p,m)))
    blocks=[d for side in sides for size in (1,2) for d in combinations(side,size)]
    M=m+8*len(blocks)
    W=1+sum(map(sum,B))
    gram=[[Q(0)]*M for _ in range(M)]
    for i,row in enumerate(B):
        for j,b in enumerate(row):
            gram[i][p+j]=gram[p+j][i]=b
    H=[]
    balance=[1]*p+[-1]*q
    for sign in (1,-1):
        h={i:sign*x for i,x in enumerate(balance)}
        for k in range(len(blocks)):
            h[m+8*k]=h[m+8*k+1]=4
        H.append(h)
    rays={}
    for k,d in enumerate(blocks):
        off=m+8*k
        for s in WORDS:
            v={off+2*i+b:1 for i,b in enumerate(s)}
            for i in d:
                if readout(s):
                    v[i]=readout(s)
            rays[k,s]=v
            for i,a in v.items():
                for j,b in v.items():
                    gram[i][j]+=W*a*b
        for i in d:
            for e,f in ((1,0),(0,1),(1,1)):
                h=defaultdict(int,{i:1})
                h[off]+=e+f
                h[off+1]+=e+f
                for j in (0,1):
                    h[off+2*j+1]-=e
                for j in (2,3):
                    h[off+2*j+1]-=f
                H.append({j:v for j,v in h.items() if v})
    E=[{i:1} for i in range(M)]+H
    # Compute E M E^t directly from the fixed, factor-independent matrix.
    target=[[sum(a*b*gram[i][j] for i,a in l.items() for j,b in r.items())
             for r in E] for l in E]
    return p,q,W,blocks,E,rays,target

def lift_test(B,atoms,require_min=True,supply=False):
    B=tuple(tuple(map(Q,row)) for row in B)
    assert sum_matrices(atoms)==B
    if require_min:
        assert rank(B)==len(atoms)  # independent certificate for the examples tested
    p,q,W,blocks,E,rays,target=model(B)
    m=p+q
    assert all(x>=0 for row in target for x in row)
    n=len(E)
    actual=[[Q(0)]*n for _ in range(n)]
    def extended(v):
        return {i:x for i,h in enumerate(E) if (x:=sum(a*v.get(j,0) for j,a in h.items()))}
    factors=[]
    core=[[Q(0)]*m for _ in range(m)]
    projected=[]
    for atom in atoms:
        a=tuple(map(sum,atom)); b=tuple(map(sum,zip(*atom))); c=sum(a)
        assert c>0 and rank(atom)==1
        z=tuple(map(Q,a+b))
        ext=extended({i:x for i,x in enumerate(z) if x})
        factors.append((Q(1,c),ext))
        recovered=tuple(tuple(Q(ext.get(i,0)*ext.get(p+j,0))/c for j in range(q)) for i in range(p))
        assert recovered==atom
        projected.append(recovered)
        for i,j in product(range(m),repeat=2):
            core[i][j]+=z[i]*z[j]/c
    for k,d in enumerate(blocks):
        if len(d)==2:
            delta=-core[d[0]][d[1]]
        else:
            i=d[0]
            side=range(p) if i<p else range(p,m)
            delta=sum(core[i][j] for j in side if j!=i)-core[i][i]
        for s in WORDS:
            weight=W+Q(delta,2)*parity(s)
            assert weight>0
            factors.append((weight,extended(rays[k,s])))
    assert len(factors)==len(atoms)+16*len(blocks)
    for weight,v in factors:
        assert min(v.values())>=0
        if supply:
            squares=list(rank_one_squares(weight))
            assert sum(x*x for x in squares)==weight
        for i,x in v.items():
            for j,y in v.items():
                actual[i][j]+=weight*x*y
    assert actual==target
    return n,len(factors),sum(actual[i][i] for i in range(n)),tuple(tuple(row) for row in core)

def shitov_cycle_test():
    # B=[1], alpha=2, root-free construction over a real quadratic field, x=y=1.
    # Independently test both orientations with rational entries for a second input.
    b=Q(100,169); alpha=isqrt(b.numerator//b.denominator)+1
    assert alpha==1
    def scalar_factor(p1,p2):
        rows=[[Q(0)]*9 for _ in range(10)]
        for which,t in enumerate([p1,p2]):
            for j in range(4):
                col=which*4+j
                rows[which][col]=t if j%2==0 else 2-t
                rows[2+4*which+j][col]=rows[2+4*which+(j+1)%4][col]=1
        rows[0][-1]=rows[1][-1]=Q(10,13)
        return rows
    u,v=scalar_factor(Q(1,13),Q(1,13)),scalar_factor(Q(25,13),Q(1,13))
    def gram(u): return [[dot(a,b) for b in u] for a in u]
    t=gram(u)
    assert sorted(zip(*u))!=sorted(zip(*v)) and t==gram(v)
    # Independent target from Shitov's definition, all 100 entries.
    target=[[Q(0)]*10 for _ in range(10)]
    target[0][0]=target[1][1]=8
    target[0][1]=target[1][0]=b
    for side in range(2):
        for j in range(4):
            i=2+4*side+j; k=2+4*side+(j+1)%4
            target[i][i]=target[side][i]=target[i][side]=2
            target[i][k]=target[k][i]=1
    assert t==target
    print("Shitov four-cycle reduction: two inequivalent rational shortest factors checked against all 100 entries.")


def dictionary_check():
    p,q,W,blocks,E,rays,target=model(((3,2),(2,3)))
    m=p+q
    def modular_rank(private):
        prime=101
        basis={}
        for v in rays.values():
            a={(i,j):x*y%prime for i,x in v.items() for j,y in v.items()
               if i<=j and (not private or j>=m)}
            while a:
                pos=min(a)
                if pos not in basis:
                    mult=pow(a[pos],-1,prime)
                    basis[pos]={k:x*mult%prime for k,x in a.items()}
                    break
                mult=a[pos]
                for k,w in basis[pos].items():
                    z=(a.get(k,0)-mult*w)%prime
                    if z:
                        a[k]=z
                    elif k in a:
                        del a[k]
        return len(basis)
    # This matrix of columns is built from raw outer products, without the moment formulas.
    assert modular_rank(True)==90
    assert modular_rank(False)==96
    print("Global auxiliary dictionary: raw upper-triangle ranks mod 101 are 90 and 96.")


def check_arithmetic_examples():
    """Verify the quintic and five-cycle examples with exact SymPy arithmetic."""
    from sympy import Matrix, Poly, QQ, Rational, eye, oo, simplify, sqrt, symbols, zeros

    x = symbols("x")
    quintic = Poly(x**5 - 4*x + 2, x, domain=QQ)
    assert quintic.is_irreducible
    assert quintic.count_roots(-oo, oo) == 3
    assert quintic.count_roots(0, 1) == 1
    assert quintic.count_roots(Rational(50849, 100000), Rational(50850, 100000)) == 1
    group, _ = quintic.galois_group()
    assert group.degree == 5 and group.order() == 120
    print("Quintic: irreducible over Q, three real roots, one in [0,1], Galois group S_5.")

    # If a Möbius transformation has matrix [[a,b],[c,d]], its finite
    # fixed points satisfy c*x**2 + (d-a)*x - b = 0 away from the pole.
    psi = Matrix([[3, -1], [1, 0]])
    a, b, c, d = list(psi**5)
    assert c != 0  # infinity is not fixed
    fixed = Poly(c*x*x + (d-a)*x - b, x, domain=QQ)
    assert fixed.monic() == Poly(x*x - 3*x + 1, x, domain=QQ)
    assert fixed.gcd(Poly(c*x + d, x, domain=QQ)).degree() == 0
    roots = [(3 - sqrt(5))/2, (3 + sqrt(5))/2]
    target = 3*eye(5)
    for i in range(5):
        j = (i+1) % 5
        target[i,j] = target[j,i] = 1
    for root in roots:
        assert root.is_positive and simplify(fixed.eval(root)) == 0
        assert simplify(3 - 1/root - root) == 0
        total = zeros(5)
        for i in range(5):
            j = (i+1) % 5
            term = zeros(5)
            term[i,i], term[j,j] = root, 1/root
            term[i,j] = term[j,i] = 1
            assert term.extract([i,j], [i,j]).det().simplify() == 0
            total += term
        assert (total - target).applyfunc(simplify) == zeros(5)
    print("Five-cycle: psi^5 has exactly two fixed points, both positive; both sums equal 3I + Adj(C_5).")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--algebra", action="store_true", help="also check the quintic and five-cycle examples using SymPy")
    args = parser.parse_args()
    if args.algebra:
        check_arithmetic_examples()
    local_test()
    dictionary_check()
    shitov_cycle_test()
    B=((3,2),(2,3))
    cores=set()
    for t in (Q(0),Q(1,4),Q(1,2),Q(2,3)):
        atoms=[outer((1,t),(3,2)),outer((0,1),(2-3*t,3-2*t))]
        n,width,_,core=lift_test(B,atoms)
        cores.add(core)
    assert len(cores)==4 and (n,width)==(78,98)
    lift_test(B,[outer((int(i==0),int(i==1)),(b*int(j==0),b*int(j==1)))
                 for i,row in enumerate(B) for j,b in enumerate(row)],require_min=False,supply=True)
    for B0,atoms in [
        (((Q(5,7),),),[outer((1,),(Q(5,7),))]),
        (((1,0,2),(0,0,0)),[outer((1,0),(1,0,2))]),
        (((1,0),(0,2),(1,1)),[outer((1,0,1),(1,0)),outer((0,2,1),(0,1))])
    ]:
        lift_test(B0,atoms,supply=True)
    print("Four-bit reduction: 7 shortest lifts, 1 nonminimal lift; all extended Gram entries and inverse cross terms agree.")
    print("One 78-by-78 matrix, four distinct rational cores, shortest width 98. Rational-square supply checks passed.")

