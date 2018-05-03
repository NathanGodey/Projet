#---------------- Création du maillage carré en 2 dimensions ------------------

"""
Marche bien concernant le calcul de codim0, codim1, codim0to1A, et codim0to1B (cf. code).
Le code calcul plein d'autres trucs potentiellement inutiles
"""

import math
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

x0 = 0; y0 = 0 # On définit l'origine du repère 2D
# Domaine de résolution
Lx = 1 # intervalle [x0,Lx] selon x
Ly = 1 # intervalle [y0,Ly] selon y

LL = min(Lx,Ly)
mylevel = 3 # Paramètre entier à modifier pour modifier le pas du maillage
            # Vaut 5 dans le code scilab
NN = 2**mylevel

Nx = NN * int(Lx/LL) # nombre de frontières selon la direction x
Ny = NN * int(Ly/LL) # nombre de frontières selon la direction y

hx = Lx/Nx # pas d'espace selon x
hy = Ly/Ny # pas d'espace selon y
volume = hx*hy # volume d'une cellule

neighbours = 4 # nombre de voisins par cellule (on met artificiellement -1 aux bords)

# Initialisation
Ncell = 0   # nombre de cellules
codim0 = np.array([[0,0,0,0,0,0,0,0]]) # chaque ligne stocke les coordonnées (x,y) du centre de chaque cellule
Nface = 0   # nombre de frontières
codim1 = np.array([[0,0,0,0]]) # chaque ligne stocke les coordonnées (x,y)
        # de deux points délimitant une frontiere 
        # (de la droite vers la gauche pour les frontières horizontales
        # du bas vers le haut pour les frontières verticales)

#Nx = Nx+1; Ny = Ny+1;
diffusion = np.zeros((Nx*Ny,Nx*Ny)) # ??? (Matrice creuse nulle)
source = np.zeros((Nx*Ny,4))
codim0to1E = []  # liste des normes des frontières (correspond avec codim1)
codim0to1NX = [] # ??? (correspond avec codim1)
codim0to1NY = [] # ??? (correspond avec codim1)
Nghost = 0 # nombre de cellules fantome?? (au bord)
codim0to1A = [] # liste des numéros des cellules en amont des frontières correspondant à codim1
codim0to1B = [] # liste des numéros des cellules en aval des frontières correspondant à codim1

fboundary = [] #???
boundary = [] #???
periodic = [] #???
antiperiodic = [] #???
for i in range(4):
        fboundary.append([])
        boundary.append([])
        periodic.append([])
        antiperiodic.append([])


# Création du maillage
# On parcours le maillage ligne par ligne, de la droite vers la gauche, et de haut en bas
ny = 0
while ny<=Ny:
    nx = 0
    while nx<=Nx:
        #(numcell = nx + ny*Nx # numéro de la cellule (on commence la numérotation à 0))
        
        # -------------- On ajoute une cellule ------------------------------
        xx = x0 + nx*hx # abscisse de la cellule
        yy = y0 + ny*hy # ordonnée de la cellule
        #print(xx,yy)
        #lt.pause(0.1)
        codim0=np.concatenate((codim0,np.array([[xx,yy,0,0,0,0,0,0]])))
        
        #diffusion[Ncell,Ncell] = 2 #?
        
        # --- On crée la frontière verticale à gauche de la cellule ---
        # Cas du bord x=0
        if nx==0:
            Nface+=1
            codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy-hy/2,xx-hx/2,yy+hy/2]])))
            plt.plot([codim1[-1,0],codim1[-1,2]],[codim1[-1,1],codim1[-1,3]],"b")
            
            ex = codim1[Nface,2]-codim1[Nface,0] # >0 (frontière horizontale), or 0
            ey = codim1[Nface,3]-codim1[Nface,1] # >0 (frontière verticale), or 0
            assert(ex>=0 and ey>=0)
            E = math.sqrt(ex*ex + ey*ey) # norme de la frontière
            NX = ey/E; #??? (ne pas confondre avec Nx!!)
            NY = ex/E; #??? (ne pas confondre avec Ny!!)
            codim0to1E.append(E)
            codim0to1NX.append(NX)
            codim0to1NY.append(NY)
            
            # On ajoute les deux côtés de la frontière
            codim0to1A.append(-(Nghost+1)) # car il n'y a pas de cellule en amont
            codim0to1B.append(Ncell) # la cellule que l'on vient de créer est en aval
            Nghost+=1
            
            fboundary[3].append(Nface) #?
            codim0[-1,3] = 4 #?
            boundary[3].append(Ncell) #?
            periodic[3].append(Ncell+(Nx-1)) #?
            antiperiodic[3].append(Nx*(Ny-(ny-1)))
            #source[Ncell,3] = 1
            codim0[-1,7] = Nface # face 4 de la cellule
        else:
            #diffusion[Ncell,Ncell-1] = -1
            codim0[Ncell,7] = Nface-2 # face 4 de la cellule
        
        # --- On crée la frontière verticale à droite de la cellule ---
        Nface+=1
        codim1=np.concatenate((codim1,np.array([[xx+hx/2,yy-hy/2,xx+hx/2,yy+hy/2]])))
        
        ex = codim1[Nface,2]-codim1[Nface,0] # >0 (frontière horizontale), or 0
        ey = codim1[Nface,3]-codim1[Nface,1] # >0 (frontière verticale), or 0
        assert(ex>=0 and ey>=0)
        E = math.sqrt(ex*ex + ey*ey) # norme de la frontière
        NX = ey/E; #??? (ne pas confondre avec Nx!!)
        NY = ex/E; #??? (ne pas confondre avec Ny!!)
        codim0to1E.append(E)
        codim0to1NX.append(NX)
        codim0to1NY.append(NY)
        
        # Cas du bord x=Nx 
        if nx==Nx:
            codim0to1A.append(Ncell) # la cellule que l'on vient de créer est en amont
            codim0to1B.append(-(Nghost+1)) # car il n'y a pas de cellule en aval
            Nghost+=1
            
            fboundary[1].append(Nface) #?
            codim0[-1,3] = 2 #?
            boundary[1].append(Ncell) #?
            periodic[1].append(Ncell-(Nx-1)) #?
            antiperiodic[1].append(Nx*(Ny-ny)+1)
            #source[Ncell,1] = 1
        else:
            codim0to1A.append(Ncell) # la cellule que l'on vient de créer est en amont
            codim0to1B.append(Ncell+1)
            
            #diffusion[Ncell,Ncell+1] = -1

        codim0[Ncell,5] = Nface # frontière 2 de la cellule

        #diffusion[Ncell,Ncell] += 2
        
        # --- On crée la frontière horizontale en-dessous de la cellule ---
        # Cas où y=0
        if ny==0:
            Nface+=1
            codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy-hy/2,xx+hx/2,yy-hy/2]])))
            
            ex = codim1[Nface,2]-codim1[Nface,0] # >0 (frontière horizontale), or 0
            ey = codim1[Nface,3]-codim1[Nface,1] # >0 (frontière verticale), or 0
            assert(ex>=0 and ey>=0)
            E = math.sqrt(ex*ex + ey*ey) # norme de la frontière
            NX = ey/E; #??? (ne pas confondre avec Nx!!)
            NY = ex/E; #??? (ne pas confondre avec Ny!!)
            codim0to1E.append(E)
            codim0to1NX.append(NX)
            codim0to1NY.append(NY)
            
            codim0to1A.append(-(Nghost+1)) # il n'y a rien en amont
            codim0to1B.append(Ncell) # la cellule que l'on vient de créer est en aval
            
            Nghost+=1
            
            fboundary[0].append(Nface) #?
            codim0[-1,2] = 1 #?
            boundary[0].append(Ncell) #?
            periodic[0].append(Ncell-(Nx-1)) #?
            antiperiodic[0].append(Nx*(Ny-ny)+1)
            source[Ncell,0] = 1

            codim0[Ncell,4] = Nface # frontière 1 de la cellule
        else:
            #diffusion[Ncell,Ncell-Nx] = -1
            codim0[Ncell,4] = Nface-(2*Nx+1)-(ny==2)*(Nx-nx) # frontière 1 de la cellule
        
        # --- On crée la frontière horizontale au-dessus de la cellule ---
        Nface+=1
        codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy+hy/2,xx+hx/2,yy+hy/2]])))
        
        ex = codim1[Nface,2]-codim1[Nface,0] # >0 (frontière horizontale), or 0
        ey = codim1[Nface,3]-codim1[Nface,1] # >0 (frontière verticale), or 0
        assert(ex>=0 and ey>=0)
        E = math.sqrt(ex*ex + ey*ey) # norme de la frontière
        NX = ey/E; #??? (ne pas confondre avec Nx!!)
        NY = ex/E; #??? (ne pas confondre avec Ny!!)
        codim0to1E.append(E)
        codim0to1NX.append(NX)
        codim0to1NY.append(NY)
        
        # Cas où y=Ny-1
        if ny==Ny:
            codim0to1A.append(Ncell) # la cellule que l'on vient de créer est en amont
            codim0to1B.append(-(Nghost+1)) # il n'y a rien en aval
            
            Nghost+=1
            
            fboundary[2].append(Nface) #?
            codim0[-1,2] = 2 #?
            boundary[2].append(Ncell) #?
            periodic[2].append(Ncell-(Nx-1)) #?
            antiperiodic[2].append(Nx*(Ny-ny)+1)
            #source[Ncell,2] = 1
        else:
            codim0to1A.append(Ncell)
            codim0to1B.append(Ncell+Nx+1)
            #diffusion[Ncell,Ncell+Nx] = -1

        codim0[Ncell,6] = Nface # frontière 3 de la cellule

        # On passe à la cellule suivante
        Ncell+=1 
        nx+=1
    ny+=1

codim0=codim0[1:,:]
codim1=codim1[1:,:]


# Visualisation et vérification du maillage:
    # En noir le maillage, en bleu le passage des frontières
for i in range(len(codim1[:,1])):
    a=[codim1[i,0],codim1[i,2]]
    b=[codim1[i,1],codim1[i,3]]
    plt.plot(a,b,"black")
    if codim0to1A[i]>=0 and codim0to1B[i]>=0:
        X=[codim0[codim0to1A[i],0],codim0[codim0to1B[i],0]]
        Y=[codim0[codim0to1A[i],1],codim0[codim0to1B[i],1]]
        plt.plot(X,Y,'b')
plt.axis('equal')
plt.show()

# Sauvegarde du maillage
np.save('codim0', codim0)
np.save('codim1', codim1)
np.save('codim0to1A', codim0to1A)
np.save('codim0to1B', codim0to1B)


