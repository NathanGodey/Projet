import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from riemann_scal import *

## Recuperation des parametres du maillage
codim0=np.load('codim0.npy')
codim1=np.load('codim1.npy')
codim0to1A=np.load('codim0to1A.npy')
codim0to1B=np.load('codim0to1B.npy')
codim0to1E=np.load('codim0to1E.npy')
codim0to1NX=np.load('codim0to1NX.npy')
codim0to1NY=np.load('codim0to1NY.npy')
x0 = 0.0; y0 = 0.0 # origine du repere 2D 
Lx = 1.0; Ly = 1.0 # domaine de resolution[x0,x0+Lx]*[y0,y0+Ly]

LL = min(Lx,Ly)
mylevel = np.load('mylevel.npy') 
NN = 2**mylevel

Nx = NN * int(Lx/LL) # nombre de faces selon la direction x
Ny = NN * int(Ly/LL) # nombre de faces selon la direction y

hx = Lx/Nx # pas d'espace selon x
hy = Ly/Ny # pas d'espace selon y
volume = hx*hy # volume d'une cellule

neighbours = 4 # nombre de voisins par cellule

## Parametres 

t = 0
Tfinal = 1 # Temps final d'une simulation 

ns = 16 # nombre de simulations
THETA=np.linspace(0,np.pi/2,ns) # differents angles de vitesse entre 0 et pi/2

normev = 0.5 # cas constant : norme de la vitesse

x0,y0=0.25,0.25 # position initiale de la gaussienne
sigma = 1/50 # largeur de la gaussienne

# Definition de la vitesse
def ux(theta):
    return normev*np.cos(theta)

def uy(theta):
    return normev*np.sin(theta)

## Initialisation

def initialisation(sigma,x0,y0):
    qini=np.zeros(codim0.shape[0])
    for i in range(len(qini)):
        x,y=codim0[i]
        qini[i]=np.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma**2))
    return qini
    
qini=initialisation(sigma,x0,y0) # Condition initiale





compteur=0
qhist=[]

## Boucle sur les differents angles de la vitesse
for theta in THETA : # Cas vitesse constante, pas uniforme sur les angles
    print(theta)
    qhist.append(qini.copy()) # Recuperation des vecteurs solutions a chaque instant de la simulation
    
    q0 = qini.copy()
    
    q1 = qini.copy()
    
    compteur+=1
    print (compteur)    
    t=0
    
    uxtheta = ux(theta)
    uytheta = uy(theta)
    
    ## Boucle temporelle
    while t<Tfinal:
        k = 1000
        flux=np.zeros(codim0.shape[0]) # Matrice des flux
        # On parcours toutes les faces
        for iface in range(codim1.shape[0]):
            i=codim0to1A[iface] # Cellule en amont
            j=codim0to1B[iface] # Cellule en aval
            if i<0 or j<0: # Conditions de bord periodiques
                continue
            # Vitesse du fluide en i et j projete selon la normale a la face
            lambdai = uxtheta*codim0to1NX[iface]+uytheta*codim0to1NY[iface] 
            lambdaj = uxtheta*codim0to1NX[iface]+uytheta*codim0to1NY[iface]
            statei = q0[i] # Concentration en i
            statej = q0[j] # Concentration en j
            # Schema de Lax-Friedriech: calcul du flux en i et en j a travers la face
            [leftf,rightf,Lambda] = RIEMANN(lambdai,statei,lambdaj,statej) 
            # Mise a jour des flux
            flux[i]+=leftf*codim0to1E[iface]
            flux[j]+=rightf*codim0to1E[iface]
            # Calcul du pas de temps associz
            k=min(k,volume/(2*(hx+hy)*Lambda))
        # Mise a jour de la concentration
        q1+=flux*(k/volume)
        qhist.append(q1.copy())
        q0 = q1.copy()
        t+=k
    # Remplissage d'une matrice compose des vecteurs solutions que l'on va rajouter a la matrice de snapshots

SNAPSHOTMATRIX=np.zeros((codim0.shape[0],len(qhist)))
for ligne in range(codim0.shape[0]):
    for colonne in range(len(qhist)):
        SNAPSHOTMATRIX[ligne,colonne]=qhist[colonne][ligne]
        
print (SNAPSHOTMATRIX.shape[0]) # Compte le nombre de snapshots
np.save('SNAPSHOTMATRIX',SNAPSHOTMATRIX) # Sauvegarde. Attention a changer le nom pour differentes simulations



## Animation 3D

for q in qhist:
    plt.clf()
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    X=np.array([codim0[i,0] for i in range(codim0.shape[0])])
    Y=np.array([codim0[i,1] for i in range(codim0.shape[0])])
    x=X.reshape(Nx+1,Ny+1)
    y=Y.reshape(Nx+1,Ny+1)
    ax.plot_wireframe(x,y,q.reshape(Nx+1,Ny+1))#,False)
    plt.pause(0.1)
