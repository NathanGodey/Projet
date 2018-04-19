import math
from scipy import sparse
import numpy as np

#-------------------------------------------------------------------------
#------------------------------ Parameters -------------------------------
#-------------------------------------------------------------------------

## vicosity=0

dimsys=1

myfilename="translatinghill"

#----- time range --------------------------------------------------------

Tinitial = 0.

Tfinal = 0.1

Nsteps = 1e4

CFLhard = True

tstpcst = False ## keep same time step as finest grid at each iteration ?

#----- grid --------------------------------------------------------------

mygeofile = "geo_square.py"

Lx = 1.
Ly = 1.

def mymask(xx,yy):
    return False

## levelmax = 7; levelmin = 4;    

levelmax = 5
levelmin = levelmax

errorcomp = True

massshow = True

deltalevel = 1 ## << to compute speed of convergence

myhist = False  ## history of levelmax: to compute order method

#-------------------------------------------------------------------------

onedx = False; onedy=False;
BC=["periodic","periodic","periodic","periodic"]

def fux(t,x,y):
    return 1.
    
def fuy(t,x,y):
    return 1.

## flot analytique : x(t,x0)=x0+t
def piedx(t,x,y):
    return x-t

## flot analytique : x(t,x0)=x0+t    
def piedy(t,x,y):
    return y-t    
    
myx0 = 0.25
myy0 = 0.25


## // x and y must have the same length
def gana(t,x,y):
    xx = piedx(t,x,y)
    yy = piedy(t,x,y)
    mywdt = 0.076
    u = math.exp(-((not onedy)*(xx-myx0)**2+(not onedx)*(yy-myy0)**2)/mywdt**2)*(((xx-myx0)**2+(yy-myy0)**2) < (0.25**2-1e-14) ) 
    return u

mymin=[0]
mymax=[1]
Tinitial = 0.
Tfinal = math.sqrt(2)/4
    
def myentropy1(state):
    return 0
    
def state2plot(i,qin):
    return qin
    
#-------------------------------------------------------------------------
    
localtest = False ## flag to perform or not an entropy test within the Riemann solver

localtest2 = False ## flag to perform or not an entropy test within the Riemann solver

globaltest = False ## flag to perform or not an entropy test after projection step

import riemann_scal

## neuman0 everywhere
def myBC(i,statein,cellcoord):
    xx = cellcoord[0]
    yy = cellcoord[1]

#-------------------------------------------------------------------------

plotevery = 1. ## << postprocessing and plot options

period = Tfinal/1 ## << convert -delay 10  img*.gif anim.gif

myplot = [True,False] ## << to plot or not : as many as dimensions dimsys+1

leg = ["c"] ## <<

sty = []
sty.append([2,-1,3]) ## <<

mygif = False ## << to animate gifs :>  convert -delay 3 img*.gif anim.gif

mypdf = False ## << to flip/flop : convert .gif -flip .jpg

mycontour = False ## << works only if ~onedx ~onedy

mycontourgif = False ## << works only if mycontour

mycontourpdf = False ## << works only if mycontour

mysleep = 10

withini = False

#-------------------------------------------------------------------------
#------------------------------ E0 Parameters -------------------------------
#-------------------------------------------------------------------------

MASS = []
L0L0 = []
L0L1 = []
L0L2 = []
L2L1 = []
L2L2 = []
for i in range(dimsys):
    MASS.append([])
    L0L0.append([])
    L0L1.append([])
    L0L2.append([])
    L2L1.append([])
    L2L2.append([])
    
levelref = levelmax+deltalevel
 
## CONVERGENCE STUDY : START BY FINEST
    
for mylevel in range(levelmax+1,levelmin,-1):
    
#-------------------------------------------------------------------------
#------------------------------ GEOMETRY ----------------------------
#-------------------------------------------------------------------------
    
    LL = min(Lx,Ly)
    NN = 2**mylevel
    if(not onedy):
        Nx = NN*math.floor(Lx/LL) ## number of cells in x direction (periodic BC truely periodic)
        hx = Lx/Nx ## cells size in x direction
    else:
        Nx = 0
        hx = 1
    if(not onedx):
        Ny = NN*math.floor(Ly/LL) ## number of cells in y direction (periodic BC truely periodic)
        hy = Ly/Ny ## cells size in x direction
    else:
        Ny = 0
        hy = 1
    
    x0 = 0
    Nx = Nx+1
    y0 = 0
    Ny = Ny+1 ## FD mesh
    
    ## 2D cartesian diffusion matrix built with mesh !! scaling
    diffusion = sparse.lil_matrix((Nx*Ny,Nx*Ny)) # (Matrice creuse nulle)
    source = sparse.lil_matrix((Nx*Ny,4)) ## OK: when there are 4 types of boundaries only
    
    """  
    *****************************************************************************
    *** 2D mesh created with the following elements (finest first) ***
    "cells" are codim0 elements (control volumes in FV method) 
    - labelled from left to right, next from bottom to top, during meshing;
    NB1: for FV discretization method, each possesses 
         1 center (bijective map similarly labelled) and 1 volume ;
    NB2: in a uniform grid all control volumes have dimensions (hx,hy)
    NB3: "cell" #(modulo(nx,Nx)+(ny-1)*Nx) bottom left vertex at ((nx-1)*hx,(ny-1)*hy)
    """
    
    Ncell = 0 ## -- number of cells = size(codim0,1)
    codim0 = [] ## -- each row stores the (x,y) coordinates of the cells centers    
                 ## i.e. v1_x,v1_y,v2_x,v2_y with [v1,v2] oriented from  
                 ## bottom to top (vertical edges) or left to right (horizontal edges)
    
    Ncellmasked = 0 ## -- number of masked cells = size(codim0masked,1)
    codim0masked = []
    
    """
    "faces" are codim1 elements (interfaces between control volumes in FV method)
    - labelled as they appear during construction of cells:
        left, right, down, top for each cells
    
    |--4--|--7--|--10-|--    3*Nx+1--|
    1  1  2  2  5  3  8 =>    Nx  3*Nx-1
    |__3__|__6__|__9__|__      3*Nx__|
    """
    
    Nface = 0 ## -- number of faces = size(codim1,1)
    codim1 = [] ## -- each row stores the (x,y) coordinates of the vertices
    
    codim0to1A = []
    codim0to1B = []
    codim0to1NX = []
    codim0to1NY = []
    codim0to1E = []

    Nghost = 0
    
    volume = hx*hy ## -- uniform here: all cells (elements) are similar rectangles
    
    neighbours = 4 ## -- uniform here (except at boundaries: we artificially put -1)
    
    diffusion = diffusion/volume
    source = source/volume
    source = source/volume ## CHECK
    
    fboundary = []
    boundary = []
    periodic = []
    antiperiodic = []
    for i in range(4):
        fboundary.append([])
        boundary.append([])
        periodic.append([])
        antiperiodic.append([])
    
    exec(mygeofile)
    
    print(Ncell,"Number of cells")
    print(Nface,"Number of faces")
    
    """
     -- to build the skyline matrix with sparse Morse storage I,J,C=k
    codim0to1 = [codim0to1;zeros(Ncell,Ncell)];
    for k=1:Nface,
      if((codim0to1I(k)>0)&(codim0to1J(k)>0)) // interior face
        codim0to1(codim0to1A(k),codim0to1B(k)) = k;  // codim0to1C(k);
        codim0to1(codim0to1B(k),codim0to1A(k)) = -k; //-codim0to1C(k);
      end
    end

//******************************************************************************
"""

    ## -- for error comp wrto analytical sol, only (fine) cell centers are needed // REMOVE

    if(not onedy):
        lx = int(2**(levelref-mylevel)) ## number of refinement in x
    else:
        lx = 1
    if(not onedx):
        ly = int(2**(levelref-mylevel)) ## number of refinement in x
    else:
        ly = 1
        
    hxref = hx/lx
    hyref = hy/ly
    XXref = []
    YYref = []
    Npatterns = 0
    for llx in range(0,lx):
        for lly in range(0,ly):
            xx0 = (-(lx-1)/2+llx)*hxref ## -(hxref/2)*(lx-1)+llx*hxref;
            yy0 = (-(lx-1)/2+llx)*hyref ## -(hyref/2)*(ly-1)+lly*hyref;
            xref=[]
            yref=[]
            for i in range(Nx):
                xref.append(x0+xx0+i*hx) ## -- periodic pattern!
                yref.append(y0+yy0+i*hx) ## -- periodic pattern!
            [YY,XX]=np.meshgrid(yref,xref)
            Npatterns = Npatterns+1
            #XXref[Npatterns] = 0
    
    
    
    
    
    
