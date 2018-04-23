#Ny, Nx, Ncell


def geosquare()
ny=1 ## we go columnwise over all rows starting from "bottom"
while ny<=Ny:
    while nx<=Nx:
        Ncell = Ncell+1 ## add a new cell labelled Ncell // (ny-1)*Nx+nx
        