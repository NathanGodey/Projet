import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

## Declaration des parametres et recuperation des donnees du maillage

nbmodes = 25 #nombre de modes pour la POD

# Domaine de resolution
Lx = 1 # intervalle [x0,Lx] selon x
Ly = 1 # intervalle [y0,Ly] selon y

LL = min(Lx,Ly)
mylevel = np.load('mylevel.npy') 
NN = 2**mylevel

Nx = NN * int(Lx/LL) # nombre de faces selon la direction x
Ny = NN * int(Ly/LL) # nombre de faces selon la direction y

hx = Lx/Nx # pas d'espace selon x
hy = Ly/Ny # pas d'espace selon y
volume = hx*hy # volume d'une cellule

codim0=np.load('codim0.npy')
codim1=np.load('codim1.npy')
codim0to1A=np.load('codim0to1A.npy')
codim0to1B=np.load('codim0to1B.npy')
codim0to1E=np.load('codim0to1E.npy')
codim0to1NX=np.load('codim0to1NX.npy')
codim0to1NY=np.load('codim0to1NY.npy')


## Recuperation de la matrice de snapshots 

M=np.load('SNAPSHOTMATRIX.npy')

## Decomposition SVD 

u,s,v=np.linalg.svd(M,False) # Thin SVD

## Trace des valeurs singulieres 
plt.plot([i for i in range(nbmodes)],[np.log(s[i]/s[0]) for i in range(nbmodes)])
plt.show()


## Trace du premier mode propre en 3D

fig=plt.figure(1)
ax=fig.gca(projection='3d')
X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
x=X.reshape(Nx+1,Ny+1)
y=Y.reshape(Nx+1,Ny+1)
ax.plot_wireframe(x,y,-u[:,0].reshape(Nx+1,Ny+1))
plt.show()

## Sauvegarde de la decomposition

np.save('u_eulerien_constant',u)
np.save('s_eulerien_constant',s)

## Recuperation et sauvegarde de la base de vecteurs propres

ur=u[:,[i for i in range(13)]]

np.save('ur_eulerien_constant',ur)


