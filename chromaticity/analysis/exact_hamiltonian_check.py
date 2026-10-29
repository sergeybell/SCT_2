# Independent exact-Hamiltonian check (fixed fields, sector bends, hard edge), delta = dp/p
import numpy as np
Br=3.47568504; L=5.3294148; h=np.radians(22.5)/L
k={'QF':1.81168058/Br,'QD':-2.01793283/Br,'QF1':1.93455213/Br,'QD1':-1.93455213/Br}
arc="OS2 QF BM QD OS1 QD BM QF SF2 QF BM QD SD QD BM QF SF1 QF BM QD SD QD BM QF SF2 QF BM QD OS1 QD BM QF OS2".split()
st="OS2 QF1 OBM QD1 OS1 QD1 OBM QF1 OS2".split()
seq=arc+st+arc+st
Ls={'BM':L,'OBM':L,'OS1':.12,'OS2':.06,'SF1':.12,'SF2':.12,'SD':.12,'QF':.21,'QD':.21,'QF1':.21,'QD1':.21}
def rhs(z,hh,K,d):
    x,px,y,py=z; ps=np.sqrt((1+d)**2-px**2-py**2)
    return np.array([(1+hh*x)*px/ps, hh*ps - hh - hh*hh*x - K*x, (1+hh*x)*py/ps, K*y])
def elem(z,e,d):
    l=Ls[e]; hh=h if e=='BM' else 0.; K=k.get(e,0.)
    if hh==0 and K==0:
        x,px,y,py=z; ps=np.sqrt((1+d)**2-px**2-py**2); return z+np.array([l*px/ps,0,l*py/ps,0])
    n=600 if hh else 60; s=l/n
    for _ in range(n):
        k1=rhs(z,hh,K,d);k2=rhs(z+s/2*k1,hh,K,d);k3=rhs(z+s/2*k2,hh,K,d);k4=rhs(z+s*k3,hh,K,d)
        z=z+s/6*(k1+2*k2+2*k3+k4)
    return z
def turn(z,d):
    for e in seq: z=elem(z,e,d)
    return z
def jac(z,d,eps=1e-7):
    f0=turn(z,d); J=np.zeros((4,4))
    for j in range(4):
        dz=np.zeros(4);dz[j]=eps; J[:,j]=(turn(z+dz,d)-turn(z-dz,d))/(2*eps)
    return f0,J
def tunes(d):
    z=np.zeros(4)
    for _ in range(4):
        f,J=jac(z,d); z=z+np.linalg.solve(np.eye(4)-J,f-z)
    f,J=jac(z,d)
    q=[]
    for i in (0,2):
        m=J[i:i+2,i:i+2]; c=(m[0,0]+m[1,1])/2; s=np.sign(m[0,1])*np.sqrt(1-c*c); q.append(np.mod(np.arctan2(s,c),2*np.pi)/2/np.pi)
    return np.array(q)
dd=1e-4
qp,q0,qm=tunes(dd),tunes(0),tunes(-dd)
print("Q0",q0); print("xi (dQ/d(dp/p))",(qp-qm)/2/dd)
