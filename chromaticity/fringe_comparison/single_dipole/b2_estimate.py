# b2_estimate.py - thin-sextupole estimate of the xi shift produced by Bell's b2 = h0 F'' theta^2 term
import numpy as np, pandas as pd
src = open('compare_single.py', encoding='utf8').read().split('cosy = {')[0]; ns = {}; exec(src, ns)
b = pd.read_csv('output/bell_bare_bm.csv'); opt = pd.read_csv('output/ring_bm_optics.csv'); k2 = pd.read_csv('output/b2_strength.csv')
out = []
for gap in (0.05, 0.025, 0.0125):
    d = b[(b.model == 'SoftEnge') & np.isclose(b.gap, gap) & (b.out <= 4)]
    m = {(r.out, (r.ex, r.epx, r.ey, r.epy, r.ed)): r.coef for r in d.itertuples()}
    R = np.array([[m.get((i+1, tuple(1 if k == j else 0 for k in range(5))), 0) for j in range(4)] for i in range(4)])
    Dc = np.array([m.get((i+1, (0, 0, 0, 0, 1)), 0) for i in range(4)])
    K = float(k2[np.isclose(k2.gap, gap)].K2L_b2.iloc[0]); sx = sy = 0.0
    for r in opt[np.isclose(opt.gap, gap)].itertuples():
        M = np.array([getattr(r, f'M{i}{j}') for i in range(1, 5) for j in range(1, 5)]).reshape(4, 4)
        D = np.array([r.D1, r.D2, r.D3, r.D4])
        for MM, DD in [(M, D), (R @ M @ np.linalg.inv(R), R @ D + Dc)]:
            bx = MM[0, 1]/np.sin(2*np.pi*r.Qx); by = MM[2, 3]/np.sin(2*np.pi*r.Qy)
            sx += K*bx*DD[0]; sy += K*by*DD[0]
    out.append(dict(gap=gap, K2L_edge=K, dxi_x_b2=sx/(4*np.pi), dxi_y_b2=-sy/(4*np.pi)))
o = pd.DataFrame(out); o.to_csv('output/b2_estimate.csv', index=False); print(o.to_string())
