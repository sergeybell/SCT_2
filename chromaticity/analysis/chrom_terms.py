import numpy as np
Br=3.47568504; L=5.3294148; h=np.radians(22.5)/L
k={'QF':1.81168058/Br,'QD':-2.01793283/Br,'QF1':1.93455213/Br,'QD1':-1.93455213/Br}
arc="OS2 QF BM QD OS1 QD BM QF SF2 QF BM QD SD QD BM QF SF1 QF BM QD SD QD BM QF SF2 QF BM QD OS1 QD BM QF OS2".split()
st="OS2 QF1 OBM QD1 OS1 QD1 OBM QF1 OS2".split(); seq=arc+st+arc+st
Ls={'BM':L,'OBM':L,'OS1':.12,'OS2':.06,'SF1':.12,'SF2':.12,'SD':.12,'QF':.21,'QD':.21,'QF1':.21,'QD1':.21}
def m2(K,l):
    if abs(K)<1e-15: return np.array([[1,l],[0,1.]])
    if K>0: s=np.sqrt(K); return np.array([[np.cos(s*l),np.sin(s*l)/s],[-s*np.sin(s*l),np.cos(s*l)]])
    s=np.sqrt(-K); return np.array([[np.cosh(s*l),np.sinh(s*l)/s],[s*np.sinh(s*l),np.cosh(s*l)]])
steps=[]
for e in seq:
    n=400 if e=='BM' else 20; hh=h if e=='BM' else 0.
    for i in range(n):
        ds=Ls[e]/n; Kx=k.get(e,0)+hh*hh; Ky=-k.get(e,0)
        A=np.eye(3); A[:2,:2]=m2(Kx,ds)
        if hh: A[0,2]=(1-np.cos(hh*ds))/hh; A[1,2]=np.sin(hh*ds)
        steps.append((e,i,n,A,m2(Ky,ds),hh,ds,Kx,Ky))
Mx=np.eye(3); My=np.eye(2)
for s in steps: Mx=s[3]@Mx; My=s[4]@My
D=np.linalg.solve(np.eye(2)-Mx[:2,:2],Mx[:2,2])
def tw(M): c=(M[0,0]+M[1,1])/2; sn=np.sign(M[0,1])*np.sqrt(1-c*c); return M[0,1]/sn,(M[0,0]-M[1,1])/(2*sn)
bx,ax=tw(Mx[:2,:2]); by,ay=tw(My)
Bx=np.array([[bx,-ax],[-ax,(1+ax*ax)/bx]]); By=np.array([[by,-ay],[-ay,(1+ay*ay)/by]]); d=np.array([D[0],D[1],1.])
T={k_:0.0 for k_ in ["quad_x","quad_y","h2bx","hDgx","hDpax","hDgy","edge_y"]}
for e,i,n,A,B,hh,ds,Kx,Ky in steps:
    if e=='BM' and i==0: T["edge_y"]+=h*By[0,0]*d[1]
    Bx2=A[:2,:2]@Bx@A[:2,:2].T; By2=B@By@B.T; d2=A@d
    av=lambda a,b:(a+b)/2
    bxm,axm,gxm=av(Bx[0,0],Bx2[0,0]),-av(Bx[0,1],Bx2[0,1]),av(Bx[1,1],Bx2[1,1])
    bym,gym=av(By[0,0],By2[0,0]),av(By[1,1],By2[1,1]); Dm,Dpm=av(d[0],d2[0]),av(d[1],d2[1])
    T["quad_x"]+=bxm*Kx*ds; T["quad_y"]+=bym*Ky*ds
    T["h2bx"]+=hh*hh*bxm*ds; T["hDgx"]+=hh*Dm*gxm*ds; T["hDpax"]+=hh*Dpm*axm*ds; T["hDgy"]+=hh*Dm*gym*ds
    Bx,By,d=Bx2,By2,d2
    if e=='BM' and i==n-1: T["edge_y"]-=h*By[0,0]*d[1]
p=4*np.pi
xi_x_lin=-T["quad_x"]/p; xi_y_lin=-T["quad_y"]/p
cx=[("A  -(1/4π)∮β_x(K1+h²)ds  = -(1/4π)∮γ_x ds  [-δp_x²/2]",xi_x_lin),
    ("C  +(1/4π)∫hDγ_x ds  [(1+hx)p_x²/2 on dispersion orbit]",T["hDgx"]/p),
    ("E  -(1/2π)∫hD'α_x ds  [cross term h D'δ x_β p_β]",-2*T["hDpax"]/p)]
cy=[("A  -(1/4π)∮β_y(-K1)ds = -(1/4π)∮γ_y ds  [-δp_y²/2]",xi_y_lin),
    ("C  +(1/4π)∫hDγ_y ds  [(1+hx)p_y²/2 on dispersion orbit]",T["hDgy"]/p),
    ("F  +(h/4π)Σ[β_y D'|in - β_y D'|out]  [hard-edge B_s kick Δp_y=∓h x' y]",T["edge_y"]/p)]
RES={"x":cx,"y":cy}
print(f"D(0)={D[0]:.6e} D'(0)={D[1]:.6e} betx0={bx:.5f} bety0={by:.5f}")
for nm,cs in (("x",cx),("y",cy)):
    s=0
    for t_,v in cs: s+=v; print(f"  xi_{nm} {t_:80s} {v:+.6f}   running {s:+.6f}")
