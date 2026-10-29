"""Analysis of the COSY OV 3 3 0 one-turn map (output/map_ov330_full.txt).
Coordinates (COSY): x, a, y, b, l, dK (dK = delta_K = dK/K)."""
import re, numpy as np, itertools, json, pathlib
P = pathlib.Path(__file__).parent/"output"
txt = (P/"map_ov330_full.txt").read_text().splitlines()
g = float(txt[0].split()[1]); C = float(txt[1].split()[1]); F = 1+1/g
maps, cur = [], None
for l in txt[2:]:
    if l.strip().startswith("I "): cur = {}; maps.append(cur); continue
    if l.strip().startswith("---") or not l.strip(): continue
    t = l.split()
    if len(t) >= 9: cur[tuple(map(int, t[3:9]))] = float(t[1])
# identity rows that COSY prints as empty
names = ["x","a","y","b","l","dK"]
def coef(i, e): return maps[i].get(tuple(e), 0.0)
def mon(e): return "".join(n*k if k<=1 else f"{n}^{k}" for n,k in zip(names,e)) if any(e) else "1"
out = []
def pr(*s):
    s = " ".join(str(v) for v in s); print(s); out.append(s)
pr(f"gamma0 = {g:.12f}, C = {C:.9f} m, 1+1/gamma0 = {F:.12f}")
# ---- 1. linear 6x6
M = np.zeros((6,6))
for i in range(6):
    for j in range(6):
        e=[0]*6; e[j]=1; M[i,j]=coef(i,e)
pr("\n## 1. Linear part R (rows: final x,a,y,b,l,dK; cols: initial)")
pr("        " + "".join(f"{n:>16s}" for n in names))
for i in range(6): pr(f"{names[i]:>6s}  " + "".join(f"{M[i,j]:16.9e}" for j in range(6)))
pr("\nkey linear terms:")
for (i,j) in [(0,5),(1,5),(4,0),(4,1),(4,5)]: pr(f"  ({names[i]}|{names[j]}) = {M[i,j]: .12e}")
pr("\n## 2. Second-order chromatic terms (i | j dK)")
for i in range(4):
    for j in range(4):
        e=[0]*6; e[j]+=1; e[5]+=1; v=coef(i,e)
        if v!=0: pr(f"  ({names[i]}|{names[j]} dK) = {v: .12e}")
for i in range(2):
    e=[0]*6; e[5]=2; pr(f"  ({names[i]}|dK dK) = {coef(i,e): .12e}")
pr("geometric 2nd-order terms used on the closed orbit:")
for i,(e) in itertools.product(range(4),[(2,0,0,0,0,0),(1,1,0,0,0,0),(0,2,0,0,0,0),(1,0,1,0,0,0),(1,0,0,1,0,0),(0,1,1,0,0,0),(0,1,0,1,0,0)]):
    v=coef(i,e)
    if v!=0: pr(f"  ({names[i]}|{mon(e)}) = {v: .12e}")
# ---- evaluate transverse map (drop l: no dependence of x,a,y,b on l w/o RF)
lterm = max(abs(v) for i in range(4) for e,v in maps[i].items() if e[4]>0) if any(e[4]>0 for i in range(4) for e in maps[i]) else 0.0
pr(f"\nmax |coef| of x,a,y,b on l: {lterm} -> l is ignorable, transverse map depends on dK only as a parameter")
def T(z, d):
    v = [z[0],z[1],z[2],z[3],0.0,d]
    return np.array([sum(c*np.prod([v[k]**e[k] for k in range(6)]) for e,c in maps[i].items()) for i in range(4)])
def J(z, d, h=1e-7):
    return np.array([(T(z+h*np.eye(4)[k],d)-T(z-h*np.eye(4)[k],d))/(2*h) for k in range(4)]).T
def co(d):
    z=np.zeros(4)
    for _ in range(50):
        dz=np.linalg.solve(np.eye(4)-J(z,d), T(z,d)-z); z=z+dz
        if np.max(abs(dz))<1e-16: break
    return z
def tune(m):
    c=(m[0,0]+m[1,1])/2; return np.arccos(c)/(2*np.pi) if m[0,1]>0 else 1-np.arccos(c)/(2*np.pi)
pr("\n## 3. Step-by-step chromaticity from the 6D map (closed orbit at fixed dK)")
D = np.linalg.solve(np.eye(2)-M[:2,:2], M[:2,5])
pr(f"linear dispersion per dK at s=0: D = {D[0]:.9e} m, D' = {D[1]:.9e}   (per dp/p: x{F:.6f} -> {D[0]*F:.9e} m)")
res={}
h=1e-4
rows=[]
for d in (-h,0.0,h):
    z=co(d); m=J(z,d); qx,qy=tune(m[:2,:2]),tune(m[2:,2:])
    rows.append((d,z,m,qx,qy))
    pr(f"dK={d:+.1e}: closed orbit x={z[0]: .9e} a={z[1]: .9e} y={z[2]:.1e}")
    pr(f"   Mx=[[{m[0,0]:.10f},{m[0,1]:.10f}],[{m[1,0]:.10f},{m[1,1]:.10f}]] Tr/2={(m[0,0]+m[1,1])/2:.12f} Qx_frac={qx:.12f}")
    pr(f"   My=[[{m[2,2]:.10f},{m[2,3]:.10f}],[{m[3,2]:.10f},{m[3,3]:.10f}]] Tr/2={(m[2,2]+m[3,3])/2:.12f} Qy_frac={qy:.12f}")
dqx=(rows[2][3]-rows[0][3])/(2*h); dqy=(rows[2][4]-rows[0][4])/(2*h)
pr(f"dQx/ddK = {dqx:.9f}, dQy/ddK = {dqy:.9f}  ->  xi_p = (1+1/g)*: xi_x = {dqx*F:.9f}, xi_y = {dqy*F:.9f}")
res["closed_orbit_FD"]=(dqx,dqy)
# analytic first derivative of trace
z1=(rows[2][1]-rows[0][1])/(2*h)
pr(f"d(closed orbit)/ddK at 0 = {z1}  (equals [D, D'] to {np.max(abs(z1[:2]-D)):.1e})")
# ---- 4. 4D map with dK as parameter: ring of 2x2 matrices M(dK) = R + dK * R1, R1 includes closed orbit feed-down
pr("\n## 4. 4D map with dK as parameter (OV 3 2 1-like extraction: drop l, dK -> parameter)")
R1=np.zeros((4,4))
for i in range(4):
    for j in range(4):
        e=[0]*6; e[j]+=1; e[5]+=1; R1[i,j]=coef(i,e)
        for k in range(2):  # feed-down from closed orbit x=D dK, a=D' dK
            e2=[0]*6; e2[j]+=1; e2[k]+=1
            fac = 2 if j==k else 1
            R1[i,j]+=fac*coef(i,e2)*D[k]
for pl,(s,nm) in enumerate([(0,"x"),(2,"y")]):
    m0=M[s:s+2,s:s+2]; m1=R1[s:s+2,s:s+2]
    c=np.trace(m0)/2; sn=np.sign(m0[0,1])*np.sqrt(1-c*c)
    dtr=np.trace(m1)
    pr(f"{nm}: dM/ddK (incl. closed-orbit feed-down) = {m1.round(9).tolist()},  dTr/ddK = {dtr:.9f}")
    q=-dtr/(4*np.pi*sn)
    pr(f"   cos mu = {c:.12f}, sin mu = {sn:.12f}, dQ/ddK = -dTr/(4 pi sin mu) = {q:.9f}, xi_p = {q*F:.9f}")
    res["param4D_"+nm]=q
    # without feed-down (what you would get forgetting the closed orbit)
    m1n=np.array([[coef(s+i,[1 if k==s+j else 0 for k in range(5)]+[1]) for j in range(2)] for i in range(2)])
    qn=-np.trace(m1n)/(4*np.pi*sn)
    pr(f"   WITHOUT closed-orbit feed-down: dTr/ddK = {np.trace(m1n):.9f} -> xi_p = {qn*F:.9f}  (feed-down term; tiny here because D(s=0)=5.9e-3 m)")
    res["nofeed_"+nm]=qn
# ---- 5. eigenvalues
pr("\n## 5. Eigenvalues of the 6x6 linear part")
ev=np.linalg.eigvals(M)
for v in ev: pr(f"  {v.real: .12f} {v.imag:+.12f}i   |lambda|={abs(v):.12f}  arg/2pi={np.angle(v)/2/np.pi: .9f}")
pr(f"longitudinal block [[ (l|l),(l|dK) ],[ (dK|l),(dK|dK) ]] = [[{M[4,4]:.6f}, {M[4,5]:.9f}], [{M[5,4]:.6f}, {M[5,5]:.6f}]]")
pr(f"(l|dK)_closed = (l|x)D+(l|a)D'+(l|dK) = {M[4,0]*D[0]+M[4,1]*D[1]+M[4,5]:.9f} m (Jordan block, eigenvalue 1 double, not diagonalisable)")
(P/"ov330_map_analysis.txt").write_text("\n".join(out)+"\n", encoding="utf-8")
json.dump({k:list(v) if isinstance(v,tuple) else v for k,v in res.items()}, open(P/"ov330_map_analysis.json","w"))
