"""rf_ts_analysis.py - can linear chromaticity be recovered from a 6D map with RF?
(1) COSY TS (normal form) on OV 3 3 0 + thin RF: dQx,y/d(r_s^2) vs V   [output/rf_ts_ov330.txt]
(2) Fixed-point method in Python: COSY OV 3 3 0 ring map (FR 0 and FR 3 dipoles) + own thin RF kick
    + constant per-turn phase slip dl (= RF frequency offset) -> off-momentum fixed point delta*(dl);
    tunes = eigen-tunes of the 6x6 Jacobian there; xi = dQ/d delta* -> compare with OV 3 2 1 as V -> 0."""
import re, json, pathlib, numpy as np
from fr_analysis import parse_maps
P = pathlib.Path(__file__).parent/"output"
res = {}
# ---------- (1) TS coefficients
t = (P/"rf_ts_ov330.txt").read_text()
ts = []
for blk in t.split("#####")[1:]:
    V = float(blk.split()[1]); body = blk.split("# TS")[1]
    first = body.split("---")[0]
    c = {tuple(map(int, e.split())): float(v) for v, e in re.findall(r"^\s+\d+\s+(\S+)\s+\d+\s+((?:\d+\s+){5}\d+)\s*$", first, re.M)}
    ts.append(dict(V_kV=V, Qx=c.get((0,)*6), dQx_drx2=c.get((2,0,0,0,0,0)), dQx_drs2=c.get((0,0,2,0,0,0)), dQx_dry2=c.get((0,0,0,0,2,0)),
                   ok=len(c) > 1))
res["TS"] = ts
# ---------- (2) fixed point with RF
def load(case):
    for b in parse_maps("map_ov330_fr.txt"):
        if b["case"][:2] == case and abs(b["case"][2]-0.05) < 1e-9: return b
def make(b):
    maps = b["maps"]; F = 1+1/b["g"]; C = b["C"]
    def ring(v):
        return np.array([sum(c*np.prod([v[k]**e[k] for k in range(6)]) for e, c in maps[i].items()) for i in range(6)])
    return ring, F, C
def rf_map(ring, C, V_kV, dl, sign):
    k = 2*np.pi/C; a = sign*V_kV*1e-3/270.0           # thin cavity, h=1: dK/K kick = a sin(k l)
    def M(v):
        w = ring(v); w[4] = w[4]+dl; w[5] = w[5]+a*np.sin(k*w[4]); return w
    return M
def jac(M, z, e=1e-7):
    return np.array([(M(z+e*np.eye(6)[i])-M(z-e*np.eye(6)[i]))/(2*e) for i in range(6)]).T
def fixed_point(M, z0):
    z = z0.copy()
    for _ in range(60):
        dz = np.linalg.solve(np.eye(6)-jac(M, z), M(z)-z); z = z+dz
        if np.max(abs(dz)) < 1e-15: break
    return z
def tunes(J):
    ev = np.linalg.eigvals(J); q = np.sort(np.abs(np.angle(ev[ev.imag > 0]))/(2*np.pi))
    return q, ev
def scan(case, Vs, h=2e-5):
    ring, F, C = make(load(case)); out = []
    for V in Vs:
        sign = None
        for s in (+1, -1):                             # choose the stable RF sign
            q, ev = tunes(jac(rf_map(ring, C, V, 0.0, s), np.zeros(6)))
            if np.all(abs(abs(ev)-1) < 1e-6) and len(q) == 3: sign = s; break
        rows = []
        for d in (-h, 0.0, h):                        # d = target fixed-point dK; dl from linear slip estimate, then solved
            M0 = rf_map(ring, C, V, 0.0, sign)
            # find dl such that fixed point has dK = d: dl = -(l(fixed transverse orbit, dK=d)) per turn
            lam = lambda dl: fixed_point(rf_map(ring, C, V, dl, sign), np.array([0,0,0,0,0,d]))
            dl = -17.8*d
            for _ in range(20):
                z = lam(dl); g = z[5]-d
                z2 = lam(dl+1e-9); dg = (z2[5]-z[5])/1e-9
                dl -= g/dg
                if abs(g) < 1e-14: break
            z = lam(dl); q, _ = tunes(jac(rf_map(ring, C, V, dl, sign), z))
            rows.append((z[5], q))
        qs = rows[1][1]; dd = rows[2][0]-rows[0][0]
        dq = (rows[2][1]-rows[0][1])/dd*F
        out.append(dict(V_kV=V, sign=sign, Q_sorted=qs.tolist(), dQ_ddelta_p_sorted=dq.tolist()))
    return out
Vs = [100, 10, 1, 0.1, 0.01, 0.001]
for case, name in [((0.0, 0.0), "FR0"), ((3.0, 0.0), "FR3dip")]:
    res["fixed_point_"+name] = scan(case, Vs)
json.dump(res, open(P/"rf_ts_analysis.json", "w"), indent=1)
for k, v in res.items():
    print("==", k)
    for r in v: print("  ", r)