#---------------- Creation du maillage carre en 2 dimensions ------------------

import math
import numpy as np
import matplotlib.pyplot as plt

x0 = 0; y0 = 0 # On definit l'origine du repere 2D (coordonnees de la premiere face)
# Domaine de resolution
Lx = 1 # intervalle [x0,Lx] selon x
Ly = 1 # intervalle [y0,Ly] selon y

LL = min(Lx,Ly)
mylevel = 8 # Parametre entier a modifier pour modifier le pas du maillage
NN = 2**mylevel

Nx = NN * int(Lx/LL) # nombre de faces selon la direction x
Ny = NN * int(Ly/LL) # nombre de faces selon la direction y

hx = Lx/Nx # pas d'espace selon x
hy = Ly/Ny # pas d'espace selon y
volume = hx*hy # volume d'une cellule

neighbours = 4 # nombre de voisins par cellule

## Initialisation
Ncell = 0   # nombre de cellules
codim0 = np.array([[0,0]]) # chaque ligne stocke les coordonnees (x,y) du centre de chaque cellule
Nface = 0   # nombre de frontieres
codim1 = np.array([[0,0,0,0]]) # chaque ligne stocke les coordonnees (x,y)
        # de deux points delimitant une face
        # (de la gauche vers la droitepour les faces horizontales
        # du bas vers le haut pour les faces verticales)

codim0to1E = []  # liste des longueurs des frontieres (correspond avec codim1)
codim0to1NX = [] # composantes x des normales des faces (correspond avec codim1)
codim0to1NY = [] # composantes y des normales des faces (correspond avec codim1)
Nghost = 0 # nombre de cellules fantome?? (au bord)
codim0to1A = [] # liste des numeros des cellules en amont des frontieres correspondant a codim1
codim0to1B = [] # liste des numeros des cellules en aval des frontieres correspondant a codim1


## Creation du maillage
# On parcours le maillage ligne par ligne, de la gauche vers la droite, et de bas en haut
ny = 0
while ny<=Ny:
    nx = 0
    while nx<=Nx:
        #(numcell = nx + ny*Nx # numero de la cellule (on commence la numerotation a 0))
        
        # -------------- On ajoute une cellule ------------------------------
        xx = x0 + nx*hx # abscisse du centre de la cellule
        yy = y0 + ny*hy # ordonnee du centre de la cellule
        
        codim0=np.concatenate((codim0,np.array([[xx,yy]])))
        
        # --- On cree la face verticale a gauche de la cellule ---
        # Cas du bord x=0
        if nx==0:
            Nface+=1
            codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy-hy/2,xx-hx/2,yy+hy/2]])))
            
            ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), ou 0
            ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), ou 0
            assert(ex>=0 and ey>=0)
            E = math.sqrt(ex*ex + ey*ey) # norme de la face
            NX = ey/E; # composante x de la normale a la face (ne pas confondre avec Nx!!)
            NY = ex/E; # composante y de la normale a la face (ne pas confondre avec Ny!!)
            codim0to1E.append(E)
            codim0to1NX.append(NX)
            codim0to1NY.append(NY)
            
            # On ajoute les deux cotes de la face
            codim0to1A.append(Ncell+Nx) # car il n'y a pas de cellule en amont
            codim0to1B.append(Ncell) # la cellule que l'on vient de creer est en aval
            Nghost+=1
            
        
        # --- On cree la face verticale a droite de la cellule ---
        Nface+=1
        codim1=np.concatenate((codim1,np.array([[xx+hx/2,yy-hy/2,xx+hx/2,yy+hy/2]])))
        
        ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), or 0
        ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), or 0
        assert(ex>=0 and ey>=0)
        E = math.sqrt(ex*ex + ey*ey) # norme de la face
        NX = ey/E; # composante x de la normale a la face (ne pas confondre avec Nx!!)
        NY = ex/E; # composante y de la normale a la face (ne pas confondre avec Ny!!)
        codim0to1E.append(E)
        codim0to1NX.append(NX)
        codim0to1NY.append(NY)
        
        # Cas du bord x=Nx 
        if nx==Nx:
            codim0to1A.append(Ncell) # la cellule que l'on vient de creer est en amont
            codim0to1B.append(-(Nghost+1)) # car il n'y a pas de cellule en aval
            Nghost+=1

        else:
            codim0to1A.append(Ncell) # la cellule que l'on vient de creer est en amont
            codim0to1B.append(Ncell+1)

        
        # --- On cree la face horizontale en-dessous de la cellule ---
        # Cas ou y=0
        if ny==0:
            Nface+=1
            codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy-hy/2,xx+hx/2,yy-hy/2]])))
            
            ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), ou 0
            ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), ou 0
            assert(ex>=0 and ey>=0)
            E = math.sqrt(ex*ex + ey*ey) # norme de la face
            NX = ey/E; # composante x de la normale a la face (ne pas confondre avec Nx!!)
            NY = ex/E; # composante y de la normale a la face (ne pas confondre avec Ny!!)
            codim0to1E.append(E)
            codim0to1NX.append(NX)
            codim0to1NY.append(NY)
            codim0to1A.append(-(Nghost+1)) # il n'y a rien en amont
            codim0to1B.append(Ncell) # la cellule que l'on vient de creer est en aval
            
            Nghost+=1
        
        # --- On cree la face horizontale au-dessus de la cellule ---
        Nface+=1
        codim1=np.concatenate((codim1,np.array([[xx-hx/2,yy+hy/2,xx+hx/2,yy+hy/2]])))
        
        ex = codim1[Nface,2]-codim1[Nface,0] # >0 (face horizontale), ou 0
        ey = codim1[Nface,3]-codim1[Nface,1] # >0 (face verticale), ou 0
        assert(ex>=0 and ey>=0)
        E = math.sqrt(ex*ex + ey*ey) # norme de la face
        NX = ey/E; # composante x de la normale a la face (ne pas confondre avec Nx!!)
        NY = ex/E; # composante y de la normale a la face (ne pas confondre avec Ny!!)
        codim0to1E.append(E)
        codim0to1NX.append(NX)
        codim0to1NY.append(NY)
        
        # Cas ou y=Ny-1
        if ny==Ny:
            codim0to1A.append(Ncell) # la cellule que l'on vient de creer est en amont
            codim0to1B.append(Ncell-Ny*(Nx+1)) # il n'y a rien en aval
            
            Nghost+=1
        else:
            codim0to1A.append(Ncell)
            codim0to1B.append(Ncell+Nx+1)

        # On passe a la cellule suivante
        Ncell+=1 
        nx+=1
    ny+=1

codim0=codim0[1:,:]
codim1=codim1[1:,:]

'''
## Visualisation et verification du maillage:
    # En noir le maillage, en bleu le passage des face
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
'''

## Sauvegarde du maillage
np.save('codim0', codim0)
np.save('codim1', codim1)
np.save('codim0to1A', codim0to1A)
np.save('codim0to1B', codim0to1B)
np.save('codim0to1E', codim0to1E)
np.save('codim0to1NX', codim0to1NX)
np.save('codim0to1NY', codim0to1NY)
np.save('mylevel', mylevel)