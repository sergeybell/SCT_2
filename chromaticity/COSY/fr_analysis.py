"""fr_analysis.py - chromaticity from OV 3 3 0 maps with fringe fields (closed-orbit method), compared with OV 3 2 1 TP.
Reads output/map_ov330_fr.txt (several CASE blocks) and output/cosy_ov321_fr.txt."""
import numpy as np, pathlib, re, json
P = pathlib.Path(__file__).parent/"output"
def parse_maps(fn):
    blocks, cur = [], None
    for l in (P/fn).read_text().splitlines():
        s = l.split()
        if not s: continue
        if s[0] == "CASE": cur = dict(case=tuple(float(v) for v in s[1:4]), maps=[]); blocks.append(cur); continue
        if s[0] == "GAMMA0": cur["g"] = float(s[1]); continue
        if s[0] == "LENGTH": cur["C"] = float(s[1]); continue
        if s[0] == "I": cur["maps"].append({}); continue
        if l.strip().startswith("---"): continue
        if len(s) >= 9: cur["maps"][-1][tuple(map(int, s[3:9]))] = float(s[1])
    return blocks
def xi_closed_orbit(b, h=1e-4):
    maps, F = b["maps"], 1+1/b["g"]
    def T(z, d):
        v = [z[0], z[1], z[2], z[3], 0.0, d]
        return np.array([sum(c*np.prod([v[k]**e[k] for k in range(6)]) for e, c in maps[i].items()) for i in range(4)])
    def J(z, d, e=1e-7):
        return np.array([(T(z+e*np.eye(4)[k], d)-T(z-e*np.eye(4)[k], d))/(2*e) for k in range(4)]).T
    def co(d):
        z = np.zeros(4)
        for _ in range(50):
            dz = np.linalg.solve(np.eye(4)-J(z, d), T(z, d)-z); z = z+dz
            if np.max(abs(dz)) < 1e-16: break
        return z
    def tune(m):
        c = (m[0,0]+m[1,1])/2; return np.arccos(c)/(2*np.pi) if m[0,1] > 0 else 1-np.arccos(c)/(2*np.pi)
    q = []
    for d in (-h, 0.0, h):
        m = J(co(d), d); q.append((tune(m[:2,:2]), tune(m[2:,2:])))
    return dict(Qx=q[1][0], Qy=q[1][1], xi_x=F*(q[2][0]-q[0][0])/(2*h), xi_y=F*(q[2][1]-q[0][1])/(2*h))
def parse_tp(fn="cosy_ov321_fr.txt"):
    t = (P/fn).read_text().split()
    out = []
    for i, w in enumerate(t):
        if w == "FD": out.append(dict(FD=float(t[i+1]), FQ=float(t[i+3]), GAP=float(t[i+5]), Qx=float(t[i+7]), Qy=float(t[i+9]), xi_x=float(t[i+11]), xi_y=float(t[i+13])))
    return out
if __name__ == "__main__":
    tp = parse_tp(); rows = []
    for b, t in zip(parse_maps("map_ov330_fr.txt"), tp):
        r = xi_closed_orbit(b)
        rows.append(dict(FR_dip=int(b["case"][0]), FR_quad=int(b["case"][1]), gap_m=b["case"][2],
                         Qx_TP=t["Qx"], Qy_TP=t["Qy"], xi_x_TP321=t["xi_x"], xi_y_TP321=t["xi_y"],
                         Qx_CO330=r["Qx"], Qy_CO330=r["Qy"], xi_x_CO330=r["xi_x"], xi_y_CO330=r["xi_y"]))
    for r in rows: print(" ".join(f"{k}={v:.9g}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))
    json.dump(rows, open(P/"fr_analysis.json", "w"), indent=1)