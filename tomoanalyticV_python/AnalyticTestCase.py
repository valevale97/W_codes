import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
from numpy.random import default_rng
from scipy.optimize import nnls

rng = default_rng(333)

def TikhonovNonNeg(W,S,v_prl,v_prp,alpha,beta):
    scale = 1/np.max(W)
    W *= scale
    
    L1_prp = (v_prp[1]-v_prp[0])**-1 * (
                np.diag(np.r_[np.zeros(v_prl.size),np.ones((v_prp.size-1)*(v_prl.size))])
                    + np.diag(-1*np.ones((v_prp.size-1)*v_prl.size),-v_prl.size))
    L1_prl = (v_prl[1]-v_prl[0])**-1 * (
                np.diag(np.tile(np.r_[-1*np.ones(v_prl.size-1),0],v_prp.size))
                + np.diag(np.tile(np.r_[1*np.ones(v_prl.size-1),0],v_prp.size)[:-1],1))
    
    indx = np.arange(v_prl.size*v_prp.size).reshape((v_prp.size,v_prl.size)).ravel(order='F')
    T = np.diag(np.zeros(v_prl.size*v_prp.size))
    for i in range(v_prl.size*v_prp.size): T[i][indx[i]] += 1
    # DFT_row = np.fft.fft(np.identity(v_prl.size*v_prp.size),norm='ortho')
    
    x_a = []
    for a in alpha:
        sL1_a = np.zeros(L1_prl.shape[1])
        sL0_a = np.zeros(W.shape[1])
        w = np.asarray([1,1,beta,0])
        # Os = 100
        # supprs = np.diag(np.r_[0,np.ones((S.size-Os)//2-1),np.zeros(Os),np.ones((S.size-Os)//2)])
        WL_a = np.concatenate([w[0]*W,w[1]*a*L1_prl,w[1]*a*L1_prp,w[2]*a*np.identity(W.shape[1])])
                                # , w[3] * supprs @ W @ DFT_row @ (np.identity(T.shape[0]) + T)
        sL_a = np.concatenate([w[0]*S,sL1_a,sL1_a,sL0_a])
                               # ,w[3] * supprs @ np.fft.fft(S,norm='ortho')
        x_a.append(nnls(WL_a,sL_a,maxiter=10**4)[0])
    
    return np.asarray(x_a)*scale

# INITIALIZING

m_p = 1.6726*10**-27
q_e = 1.6021917*10**-19
m_D = 2*m_p

T_prl = 20 # KeV
T_prp = 20 # KeV
n_e = 10**19 # m**-3
v_drift = 5*10**5 # m/s

v_th_prl = (2*T_prl*q_e*10**3 / m_D)**0.5
v_th_prp = (2*T_prp*q_e*10**3 / m_D)**0.5

n,m = 40,20
v_prl = np.linspace(-4*10**6,4*10**6,n)
v_prp = np.linspace(10**4,4*10**6,m)
v = np.meshgrid(v_prl,v_prp)

Fv_3DbiMax = n_e/(np.pi**1.5*v_th_prl*v_th_prp**2) * (np.exp(-((v[0]-v_drift)/v_th_prl)**2)*
                                                      np.exp(-((v[1])/v_th_prp)**2) )
Fv_2DbiMax = Fv_3DbiMax * 2*np.pi*v[1]

phi = np.array([10,20,40,70,85])*np.pi/180 #10,20,40,70,85

# COMPUTING (ANALYTICAL) WEIGHT FUNCTION

u = np.arange(-1,1,0.02) * 5*10**6
du = u[1]-u[0]
dv_prl,dv_prp = v_prl[1]-v_prl[0],v_prp[1]-v_prp[0]

W_synth = np.zeros((np.size(phi)*np.size(u),np.size(v[0])))

for i in range(len(phi)):
    p_i = phi[i]
    for j in range(len(u)):
        u_j = u[j]
        g1 = np.arccos((u_j-0.5*du-v[0]*np.cos(p_i)) / (v[1]*np.sin(p_i)), dtype=complex)
        g2 = np.arccos((u_j+0.5*du-v[0]*np.cos(p_i)) / (v[1]*np.sin(p_i)), dtype=complex)
        W_synth[len(u)*i+j] += np.real(g1-g2).ravel() / (np.pi*du) * (dv_prl*dv_prp)

W_synth = np.delete(W_synth,(np.sum(W_synth,axis=1)<np.finfo(float).eps),axis=0)

# ADDING NOISE AND PERFORMING INVERSION

S_synth = np.einsum('ij,j->i',W_synth,Fv_2DbiMax.ravel())

S = S_synth + 0.01*np.mean(S_synth**0.5)*(rng.standard_normal(np.size(S_synth))*
                                          np.max([S_synth**0.5,np.ones(np.size(S_synth))*1],axis=0))
e = 0.01*np.mean(np.abs(S)**0.5)*np.max([np.abs(S)**0.5,np.ones(np.size(S))*1],axis=0)

W = np.einsum('ij,i->ij',W_synth,1/e)
S = np.einsum('i,i->i',S,1/e)

# PERFORMING INVERSION

lmbd = np.asarray([10**1,10**2,10**3,10**4])#**0.5
bt = 2*10**-15

f_lambda = TikhonovNonNeg(W,S,v_prl,v_prp,lmbd,bt)


# PLOTTING RESULTS

plt.rcParams['axes.grid'] = False
plt.rcParams.update({'figure.autolayout': False})
matplotlib.rcParams.update({'font.size': 10})

k = matplotlib.cm.get_cmap('plasma_r').copy()
k.set_under('w')
normalize = matplotlib.colors.Normalize(vmin=10**2, vmax=0.75*Fv_2DbiMax.max())
kmap = matplotlib.cm.ScalarMappable(norm=normalize, cmap=k)

fig,mats = plt.subplot_mosaic("AB;CD;EF", figsize=(9,9), sharex=True,sharey=True)
fig.subplots_adjust(0.1,0.075,0.85,0.95,hspace=0.15,wspace=0.20)

mats["A"].contourf(v[0],v[1],Fv_2DbiMax,norm=normalize,cmap=k,extend='min')
mats["A"].set_title("True F")
mats["B"].contourf(v[0],v[1],np.sum(W_synth,axis=0).reshape(v[0].shape),
                   norm=None,levels=30)
mats["B"].set_title("Synthetic W")
mats["A"].set_ylabel("$v_\perp\:[m^1s^{-1}]$")

for mx in range(0,4):
    mats[list(mats)[mx+2]].contourf(v[0],v[1],f_lambda[mx].reshape(v[0].shape),
                                   norm=normalize,levels=30,cmap=k,extend='min')
    mats[list(mats)[mx+2]].set_title(f"$\lambda=\:{lmbd[mx]:.1e}$")
mats["C"].set_ylabel("$v_\perp\:[m^1s^{-1}]$")
mats["E"].set_ylabel("$v_\perp\:[m^1s^{-1}]$")
mats["E"].set_xlabel("$v_\parallel\:[m^1s^{-1}]$")
mats["F"].set_xlabel("$v_\parallel\:[m^1s^{-1}]$")

fig.colorbar(kmap,cax=fig.add_axes([0.9, 0.2, 0.025, 0.6]),label='f [$s^1m^{-4}$]')
plt.show()

#%%

# COMPARE WITH MATLAB OUTPUT

def rereshape(array,v):
    for j in range(array.shape[0]):
        array[j] = array[j,:].reshape(v.shape,order='F').ravel()
    
    return array

f_lambda_MAT = rereshape(np.genfromtxt("xalpha.dat",usecols=(0,1,2,3),unpack=True),v[0])
W_MAT = rereshape(np.genfromtxt("WS_noised.dat")[:,:-1],v[0])
S_MAT = np.genfromtxt("WS_noised.dat")[:,-1]

lmbd = np.asarray([10**1,10**2,10**3,10**4])#**0.5
bt = 2*10**-6

f_lambda_PY = TikhonovNonNeg(W_MAT,S_MAT,v_prl,v_prp,lmbd,bt)

h = matplotlib.cm.get_cmap('bwr').copy()
normalize = matplotlib.colors.Normalize(vmin=-10**-8, vmax=10**-8)
hmap = matplotlib.cm.ScalarMappable(norm=normalize,cmap=h)
nermalize = matplotlib.colors.Normalize(vmin=10**2, vmax=0.75*Fv_2DbiMax.max())
kmap = matplotlib.cm.ScalarMappable(norm=nermalize, cmap=k)

fig,mats = plt.subplot_mosaic("AB;CD;EF;GH", figsize=(9,9), sharex=True,sharey=True)
fig.subplots_adjust(0.1,0.075,0.8,0.95,hspace=0.2,wspace=0.15)

for mx in range(0,4):
    mats[list(mats)[2*mx]].contourf(v[0],v[1],f_lambda_PY[mx].reshape(v[0].shape),
                                   norm=nermalize,levels=30,cmap=k,extend='min')
    mats[list(mats)[2*mx]].set_title(f"$\lambda=\:{lmbd[mx]:.1e}$")
    mats[list(mats)[2*mx+1]].contourf(v[0],v[1],(f_lambda_PY[mx]-f_lambda_MAT[mx]).reshape(v[0].shape),
                                   norm=normalize,levels=50,cmap=h)
    mats[list(mats)[2*mx+1]].set_title("vs MATLAB")
mats["A"].set_ylabel("$v_\perp\:[m^1s^{-1}]$")
mats["C"].set_ylabel("$v_\perp\:[m^1s^{-1}]$")
mats["E"].set_ylabel("$v_\perp\:[m^1s^{-1}]$")
mats["G"].set_ylabel("$v_\perp\:[m^1s^{-1}]$")
mats["G"].set_xlabel("$v_\parallel\:[m^1s^{-1}]$")
mats["H"].set_xlabel("$v_\parallel\:[m^1s^{-1}]$")

fig.colorbar(kmap,cax=fig.add_axes([0.85, 0.55, 0.025, 0.4]),label='f [$s^1m^{-4}$]')
cb = fig.colorbar(hmap,cax=fig.add_axes([0.85, 0.05, 0.025, 0.4]),label='f [$s^1m^{-4}$]')
# cb.formatter.set_powerlimits((0, 0))
plt.show()





