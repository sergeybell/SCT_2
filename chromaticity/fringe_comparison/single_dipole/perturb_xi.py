# perturb_xi.py - attribute ring xi difference (Bell SoftEdge Enge vs COSY FR3) to individual 2nd-order
# coefficients of the bare-BM map, by first-order perturbation of the one-turn trace on the closed orbit.
import numpy as np, pandas as pd
import importlib.util, sys
spec = importlib.util.spec_from_file_location('cs', 'compare_single.py')
src = open('compare_single.py', encoding='utf8').read().split('cosy = {')[0]
ns = {}; exec(src, ns); read_cosy, to_delta = ns['read_cosy'], ns['to_delta']
cos = {0.0: to_delta(read_cosy('output/mapbm_1.txt')), 0.05: to_delta(read_cosy('output/mapbm_2.txt')),
       0.025: to_delta(read_cosy('output/mapbm_3.txt')), 0.0125: to_delta(read_cosy('output/mapbm_4.txt'))}
b = pd.read_csv('output/bell_bare_bm.csv')
def bell(model, gap):
    d = b[(b.model == model) & np.isclose(b.gap, gap) & (b.out <= 4)]
    return {(r.out, (r.ex, r.epx, r.ey, r.epy, r.ed)): r.coef for r in d.itertuples()}
bel = {0.0: bell('NoEdge', 0.0), **{g: bell('SoftEnge', g) for g in (0.05, 0.025, 0.0125)}}
opt = pd.read_csv('output/ring_bm_optics.csv')
V = ['x', 'a', 'y', 'b', 'd']
def lin(m):
    R = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            e = [0]*5; e[j] = 1; R[i, j] = m.get((i+1, tuple(e)), 0)
    return R
def jprime_terms(m, w):
    """list of (i,j,term,contribution) to dJ_ij/ddelta on closed orbit z = w*delta"""
    out = []
    for (i, e), c in m.items():
        if sum(e) != 2: continue
        for j in range(4):
            if e[j] == 0: continue
            e2 = list(e); e2[j] -= 1          # derivative wrt z_j
            fac = e[j]
            k = e2.index(1)                    # remaining variable
            out.append((i-1, j, (i, e), c*fac*w[k]))
    return out
def lab(i, e): return '(' + ['x', 'a', 'y', 'b'][i-1] + '|' + ''.join(v*k if k == 1 else v+'^2' for v, k in zip(V, e) if k) + ')'
rows = []; summary = []
for gap in (0.05, 0.025, 0.0125):
    o = opt[np.isclose(opt.gap, gap)]
    for refname, A, B in [('BellSoft-COSYFR3', bel[gap], cos[gap]), ('BellSoft-BellNoEdge', bel[gap], bel[0.0]), ('COSYFR3-COSYFR0', cos[gap], cos[0.0])]:
        Re = lin(bel[gap]); tot = {'x': 0.0, 'y': 0.0}; per = {}
        for r in o.itertuples():
            M = np.array([getattr(r, f'M{i}{j}') for i in range(1, 5) for j in range(1, 5)]).reshape(4, 4)
            w = [r.D1, r.D2, r.D3, r.D4, 1.0]
            Rr = M @ np.linalg.inv(Re)
            Q = {'x': r.Qx, 'y': r.Qy}
            for (src_m, sgn) in [(A, 1), (B, -1)]:
                for i, j, t, c in jprime_terms(src_m, w):
                    pl = 'x' if i < 2 else 'y'
                    if (j < 2) != (i < 2): continue
                    # d tr(block) = sum_ij Rr_ji dJ_ij
                    contrib = -sgn*Rr[j, i]*c/(4*np.pi*np.sin(2*np.pi*Q[pl]))
                    tot[pl] += contrib; per[(pl, t)] = per.get((pl, t), 0) + contrib
        summary.append(dict(gap=gap, pair=refname, dxi_x_pert=tot['x'], dxi_y_pert=tot['y']))
        for (pl, t), v in per.items():
            rows.append(dict(gap=gap, pair=refname, plane=pl, term=lab(*t), contrib=v))
S = pd.DataFrame(summary); T = pd.DataFrame(rows)
S.to_csv('output/perturb_summary.csv', index=False); T.to_csv('output/perturb_terms.csv', index=False)
print(S.to_string())
t = T[(T.gap == 0.05) & (T.pair == 'BellSoft-COSYFR3')].copy(); t['abs'] = t.contrib.abs()
print(t.sort_values('abs', ascending=False).head(14).to_string())
