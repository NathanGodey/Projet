ny = 1;             // -- we go columnwise over all rows starting from "bottom"
while ny<=Ny,
  nx = 1;
  while nx<=Nx,
  
  
    Ncell = Ncell+1;                   // add a new cell labelled Ncell // (ny-1)*Nx+nx
    xx = x0+(nx-1)*hx;
    yy = y0+(ny-1)*hy;
    codim0 = [codim0;[xx,yy,zeros(1,6)]];// add a new cell center to list
    
    diffusion(Ncell,Ncell) = 0;
    
  if (~onedy)  

    diffusion(Ncell,Ncell) = diffusion(Ncell,Ncell) + 2;
    
    if (nx==1) then

      Nface = Nface+1; // add a new  VERTICAL face labelled Nface
      codim1 = [codim1;[xx-hx/2,yy-hy/2,xx-hx/2,yy+hy/2] ];
      // -- tangent to interf of length E (scalar coeff 2|c_ij| of flux) --
      ex = codim1(Nface,3)-codim1(Nface,1);     // >0 (horizontal edge), or 0
      ey = codim1(Nface,4)-codim1(Nface,2);     // >0 (vertical edge), or 0
      if((ex<0)|(ey<0)) disp("Warning"); end    // sign:(nx,ny) from I=A to J=B
      E = sqrt(ex^2+ey^2);
      NX = ey/E;
      NY = ex/E;
      codim0to1E = [codim0to1E,E];
      codim0to1NX = [codim0to1NX,NX];
      codim0to1NY = [codim0to1NY,NY];
      // -- normal to interf i->j = A->B (flux sign: left/right Riemann)
      
      codim0to1A = [codim0to1A,-(Nghost+1)];
      codim0to1B = [codim0to1B,Ncell];

      Nghost = Nghost+1;
      
      fboundary(4) = [fboundary(4),Nface];
      
      codim0(Ncell,4) = 4;                             // boundary(4) flag
      boundary(4) = [boundary(4),Ncell];
      periodic(4) = [periodic(4),Ncell+(Nx-1)];            // periodic cell matching boundary(4)
      antiperiodic(4) = [antiperiodic(4),Nx*(Ny-(ny-1))];  // antiperiodic boundary(4) match
      
      source(Ncell,4) = 1;
      
      codim0(Ncell,8) = Nface;                         // face 4 of the cell

    else 
      
      diffusion(Ncell,Ncell-1) = -1;
    
      codim0(Ncell,8) = Nface-2; // face 4 of the cell ************************************
      
    end

    
    Nface = Nface+1; // add a new VERTICAL face labelled Nface
    codim1 = [codim1;[xx+hx/2,yy-hy/2,xx+hx/2,yy+hy/2] ];
    // -- tangent to interf of length E (scalar coeff 2|c_ij| of flux) --
    ex = codim1(Nface,3)-codim1(Nface,1);     // >0 (horizontal edge), or 0
    ey = codim1(Nface,4)-codim1(Nface,2);     // >0 (vertical edge), or 0
    if((ex<0)|(ey<0)) disp("Warning"); end    // sign:(nx,ny) from I=A to J=B
    E = sqrt(ex^2+ey^2);
    NX = ey/E;
    NY = ex/E;
    codim0to1E = [codim0to1E,E];
    codim0to1NX = [codim0to1NX,NX];
    codim0to1NY = [codim0to1NY,NY];
    // -- normal to interf i->j = A->B (flux sign: left/right Riemann)      
    
    if (nx==Nx) then
    
      codim0to1A = [codim0to1A,Ncell];
      codim0to1B = [codim0to1B,-(Nghost+1)];

      Nghost = Nghost+1;
      
      fboundary(2) = [fboundary(2),Nface];
      
      codim0(Ncell,4) = 2;                             // boundary(2) flag
      boundary(2) = [boundary(2),Ncell];
      periodic(2) = [periodic(2),Ncell-(Nx-1)];            // periodic cells matching boundary(2)
      antiperiodic(2) = [antiperiodic(2),Nx*(Ny-ny)+1];  // antiperiodic boundary(2) match
      
      source(Ncell,2) = 1;
      
    else  

      codim0to1A = [codim0to1A,Ncell];
      codim0to1B = [codim0to1B,Ncell+1];
      
      diffusion(Ncell,Ncell+1) = -1;
    
    end

    codim0(Ncell,6) = Nface; // face 2 of the cell ****************************************

  end
    
    // right (vertical) cell face above //////// top (horizontal) cell face below        
  
  if (~onedx)  

    diffusion(Ncell,Ncell) = diffusion(Ncell,Ncell) + 2;

    
    if (ny==1) then

      Nface = Nface+1; // add a new HORIZONTAL face labelled Nface
      codim1 = [codim1;[xx-hx/2,yy-hy/2,xx+hx/2,yy-hy/2]];
      // -- tangent to interf of length E (scalar coeff 2|c_ij| of flux) --
      ex = codim1(Nface,3)-codim1(Nface,1);   // >0 (horizontal edge), or 0
      ey = codim1(Nface,4)-codim1(Nface,2);   // >0 (vertical edge), or 0
      if((ex<0)|(ey<0)) disp("Warning"); end  // sign : (nx,ny) goes from I to J
      E = sqrt(ex^2+ey^2);      
      NX = ey/E;
      NY = ex/E;    
      codim0to1NX = [codim0to1NX,NX];
      codim0to1NY = [codim0to1NY,NY];
      codim0to1E = [codim0to1E,E];
      // -- normal to interface i->j = 1->2 (flux sign: left/right Riemann)
      
      codim0to1A = [codim0to1A,-(Nghost+1)];
      codim0to1B = [codim0to1B,Ncell];

      Nghost = Nghost+1;
      
      fboundary(1) = [fboundary(1),Nface];
      
      codim0(Ncell,3) = 1;                          // boundary flag
      boundary(1) = [boundary(1),Ncell];
      periodic(1) = [periodic(1),nx+Nx*(Ny-1)];        // periodic cells matching boundary(1)
      antiperiodic(1) = [antiperiodic(1),Nx*Ny-nx+1];  // antiperiodic boundary(1) match
      
      codim0(Ncell,5) = Nface;                     // face 1 of the cell
      
      source(Ncell,1) = 1;
      
    else
      
      diffusion(Ncell,Ncell-Nx) = -1;
    
      codim0(Ncell,5) = Nface-(2*Nx+1)-(ny==2)*(Nx-nx); // face 1 of cell ********

    end  

    
    Nface = Nface+1; // add a new HORIZONTAL face labelled Nface
    codim1 = [codim1;[xx-hx/2,yy+hy/2,xx+hx/2,yy+hy/2]];
    // -- tangent to interf of length E (scalar coeff 2|c_ij| of flux) --
    ex = codim1(Nface,3)-codim1(Nface,1);   // >0 (horizontal edge), or 0
    ey = codim1(Nface,4)-codim1(Nface,2);   // >0 (vertical edge), or 0
    if((ex<0)|(ey<0)) disp("Warning"); end  // sign : (nx,ny) goes from I to J
    E = sqrt(ex^2+ey^2);      
    NX = ey/E;
    NY = ex/E;    
    codim0to1NX = [codim0to1NX,NX];
    codim0to1NY = [codim0to1NY,NY];
    codim0to1E = [codim0to1E,E];
    // -- normal to interface i->j = 1->2 (flux sign: left/right Riemann)
    
    if (ny==Ny) then
    
      codim0to1A = [codim0to1A,Ncell];
      codim0to1B = [codim0to1B,-(Nghost+1)];
      
      Nghost = Nghost+1;
      
      fboundary(3) = [fboundary(3),Nface];
      
      codim0(Ncell,3) = 3;                               // boundary(3) flag
      boundary(3) = [boundary(3),Ncell];
      periodic(3) = [periodic(3),Ncell-Nx*(Ny-1)];       // periodic cells matching boundary(2)
      antiperiodic(3) = [antiperiodic(3),Nx*(Ny-ny)+1];  // antiperiodic boundary(2) match
      
      source(Ncell,3) = 1;
      
    else
    
      codim0to1A = [codim0to1A,Ncell];
      codim0to1B = [codim0to1B,Ncell+Nx];    
      
      diffusion(Ncell,Ncell+Nx) = -1;

    end

    codim0(Ncell,7) = Nface; // face 3 of the cell **************************

  end  
    
    nx = nx+1;
  end ///////////////////////////////////////// -- enf of while nx<=Nx -- // 
  ny = ny+1;
end /////////////////////////////////////////// -- enf of while ny<=Ny -- // 
