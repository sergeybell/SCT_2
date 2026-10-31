# compare_single.py - COSY vs Bell maps of a lone BM; writes output/single_bm_compare.csv
import numpy as np, pandas as pd, itertools
from math import sqrt
M0 = 1875.61294257; g = 1 + 270/M0; bg = sqrt(g*g-1)
def dK_series(N=3):
    a, b, c = 1+bg**2, 2*bg**2, bg**2; P = [a, b, c]+[0]*N; s = [sqrt(a)]+[0]*N
    for n in range(1, N+1): s[n] = (P[n]-sum(s[i]*s[n-i] for i in range(1, n)))/(2*s[0])
    f = np.array(s)/(g-1); f[0] -= 1/(g-1)+1; return f
fK = dK_series()
def read_cosy(fn):
    m = {}
    for L in open(fn):
        if len(L) < 75 or not L[71:].strip().isdigit(): continue
        ex = L[71:].strip()
        vals = [float(L[1+i*14:15+i*14].replace('D','E')) for i in range(5)]
        e = tuple(int(ch) for ch in ex[:4]) + (int(ex[6]),)   # x a y b dK
        for i, v in enumerate(vals[:4]):
            if v: m[(i+1, e)] = v
    return m
def to_delta(m):   # substitute dK = f1 d + f2 d^2 + f3 d^3, keep order <= 3
    out = {}
    for (i, e), v in m.items():
        k = e[4]; poly = np.array([1.0]); 
        for _ in range(k): poly = np.convolve(poly, fK)[:4]
        for n, c in enumerate(poly):
            if c == 0 or sum(e[:4])+n > 3 or (k > 0 and n == 0): continue
            key = (i, e[:4]+(n,)); out[key] = out.get(key, 0)+v*c
    return out
cosy = {'FR0': to_delta(read_cosy('output/map_1.txt'))}
for j, G in [(2, 0.05), (3, 0.025), (4, 0.0125)]: cosy[f'FR3_{G}'] = to_delta(read_cosy(f'output/map_{j}.txt'))
b = pd.read_csv('output/bell_single_bm.csv')
def bell(model, gap=0.0):
    d = b[(b.model == model) & np.isclose(b.gap, gap) & (b.out <= 4)]
    return {(r.out, (r.ex, r.epx, r.ey, r.epy, r.ed)): r.coef for r in d.itertuples()}
B = {'NoEdge': bell('NoEdge'), 'HardEdge3': bell('HardEdge3')}
for G in (0.05, 0.025, 0.0125): B[f'Soft_{G}'] = bell('SoftEnge', G); B[f'FINT_{G}'] = bell('HardEdge2_FINT', G)
names = ['x', 'px', 'y', 'py']
def lab(i, e):
    v = ''.join(n*k if k == 1 else f'{n}^{k}' for n, k in zip(['x', 'a', 'y', 'b', 'd'], e) if k)
    return f'({names[i-1]}|{v})'
keys = sorted(set().union(*[set(m) for m in list(cosy.values())+list(B.values())]), key=lambda k: (sum(k[1]), k))
rows = []
for k in keys:
    r = dict(term=lab(*k), order=sum(k[1]), out=k[0], cosy_FR0=cosy['FR0'].get(k, 0), bell_NoEdge=B['NoEdge'].get(k, 0), bell_Hard3=B['HardEdge3'].get(k, 0))
    for G in (0.05, 0.025, 0.0125):
        r[f'cosy_FR3_{G}'] = cosy[f'FR3_{G}'].get(k, 0); r[f'bell_Soft_{G}'] = B[f'Soft_{G}'].get(k, 0); r[f'bell_FINT_{G}'] = B[f'FINT_{G}'].get(k, 0)
        r[f'diff_{G}'] = r[f'bell_Soft_{G}'] - r[f'cosy_FR3_{G}']
    rows.append(r)
df = pd.DataFrame(rows); df.to_csv('output/single_bm_compare.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
d = df.copy(); d['d0'] = d.bell_NoEdge - d.cosy_FR0
d['r'] = d['diff_0.05']/d['diff_0.025'].replace(0, np.nan)
print(d[(d['diff_0.05'].abs() > 1e-7) | (d.d0.abs() > 1e-7)][['term', 'cosy_FR0', 'd0', 'cosy_FR3_0.05', 'bell_Soft_0.05', 'diff_0.05', 'diff_0.025', 'r']].to_string())
