clear;

// --------------------------------------------------------------------
// --------------------------- PARAMETERS -----------------------------
// --------------------------------------------------------------------

// viscosity = 0;

dimsys = 1;

myfilename = "translatinghill";

// --  time range -----------------------------------------------------

Tinitial = 0.;

Tfinal = .1;

Nsteps = 10^4;

CFLhard = %T; //%F

tstpcst = %F; //%F // -- keep same time step as finest grid at each iteration ?

// -- grid ----------------------------------------------------------------

mygeofile = "geo_square.sci";

Lx = 1.; Ly = 1.;

deff('[myres]=mymask(xx,yy)','myres=%F');

//levelmax = 7; levelmin = 4;

levelmax = 5; levelmin = levelmax;

errorcomp = %T;

massshow = %T;

deltalevel = 1; // << to compute speed of convergence

myhist = %F;  // -- history of levelmax: to compute order method

//*************************************************************************    

onedx = %F; onedy = %F; BC = ['periodic','periodic','periodic','periodic'];
deff('[u]=fux(t,x,y)','u=1.');
deff('[u]=fuy(t,x,y)','u=1.');
deff('[u]=piedx(t,x,y)','u=x-t'); // flot analytique : x(t,x0)=x0+t
deff('[u]=piedy(t,x,y)','u=y-t'); // flot analytique : x(t,x0)=x0+t
// A VALIDER
// deff('[u]=fux(t,x,y)','u=cos(t)');
// deff('[u]=fuy(t,x,y)','u=sin(t)');
// deff('[u]=piedx(t,x,y)','u=x-sin(t)'); // flot analytique : x(t,x0)=x0+sin(t)
// deff('[u]=piedy(t,x,y)','u=y+cos(t)-1');
// deff('[u]=masource(t,x,y,q)','u=0');
myx0 = .25;
myy0 = .25;
function u=qana(t,x,y) // x and y must have the same length
  xx = piedx(t,x,y);
  yy = piedy(t,x,y);
  mywdt = .076;
//  u = exp(-((~onedy)*(xx-myx0).^2+(~onedx)*(yy-myy0).^2)/mywdt^2);
  u = exp(-((~onedy)*(xx-myx0).^2+(~onedx)*(yy-myy0).^2)/mywdt^2)...
        .*(((xx-myx0).^2+(yy-myy0).^2)<.25^2-1e-14);
endfunction
mymin = [0];
mymax = [1];
Tinitial = 0.; // <<
Tfinal = sqrt(2)/4; // <<

deff('[myres]=myentropy1(state)','myres=0');

function qout=state2plot(i,qin)
  qout=qin
endfunction

//*************************************************************************    

localtest = %F; // ** flag to perform or not an entropy test within the Riemann solver

localtest2 = %F; // ** flag to perform or not an entropy test within the Riemann solver

globaltest = %F; // ** flag to perform or not an entropy test after projection step

exec("riemann_scal.sci");

function stateout = myBC(i,statein,cellcoord) // neuman0 everywhere
  xx = cellcoord(1);
  yy = cellcoord(2);
endfunction 

//*************************************************************************    

plotevery = 1; // << postprocessing and plot options

period = Tfinal/1; // << convert -delay 10  img*.gif anim.gif
 
myplot = [%T,%F]; // << to plot or not : as many as dimensions dimsys+1

leg = ["c"]; // <<

sty = list();
sty(1) = [2,-1,3]; // <<

mygif = %F;  // << to animate gifs :>  convert -delay 3 img*.gif anim.gif

mypdf = %F;  // << to flip/flop : convert .gif -flip .jpg

mycontour = %F; // << works only if ~onedx ~onedy

mycontourgif = %F; // << works only if mycontour

mycontourpdf = %F; // << works only if mycontour

mysleep = 10;

withini = %F;

// --------------------------------------------------------------------
// -------------------------- EO PARAMETERS ---------------------------
// --------------------------------------------------------------------

MASS = list();
L0L0 = list();
L0L1 = list();
L0L2 = list();
L2L1 = list();
L2L2 = list();
for i=1:dimsys
  MASS(i) = [];
  L0L0(i) = [];
  L0L1(i) = [];
  L0L2(i) = [];
  L2L1(i) = [];
  L2L2(i) = [];
end

levelref = levelmax+deltalevel;

for mylevel=levelmax:-1:levelmin //// CONVERGENCE STUDY : START BY FINEST //////

//   mylevel=levelmin

// -----------------------------------------------------------------------------
// --------------------------------- GEOMETRY ----------------------------------
// -----------------------------------------------------------------------------

LL = min(Lx,Ly);
NN = (2^mylevel);

if(~onedy)
  Nx = NN*floor(Lx/LL); // number of cells in x direction (periodic BC truely periodic)
  hx = Lx/Nx; // cells size in x direction
else
  Nx = 0;
  hx = 1;
end
if(~onedx)
  Ny = NN*floor(Ly/LL); // number of cells in y direction (periodic BC truely periodic)
  hy = Ly/Ny; // cells size in y direction
else
  Ny = 0;
  hy = 1;
end
//
x0 = 0; Nx = Nx+1; y0 = 0; Ny = Ny+1; // FD mesh
//
// x0 = hx/2; y0 = hy/2; // FV mesh

diffusion = spzeros(Nx*Ny,Nx*Ny); // 2D cartesian diffusion matrix built with mesh !! scaling
source = spzeros(Nx*Ny,4); // OK: when there are 4 types of boundaries only

// *****************************************************************************
// *** 2D mesh created with the following elements (finest first) ***
// "cells" are codim0 elements (control volumes in FV method) 
// - labelled from left to right, next from bottom to top, during meshing;
// NB1: for FV discretization method, each possesses 
//      1 center (bijective map similarly labelled) and 1 volume ;
// NB2: in a uniform grid all control volumes have dimensions (hx,hy)
// NB3: "cell" #(modulo(nx,Nx)+(ny-1)*Nx) bottom left vertex at ((nx-1)*hx,(ny-1)*hy)
//
Ncell = 0;          // -- number of cells = size(codim0,1)
codim0 = [];        // -- each row stores the (x,y) coordinates of the cells centers
//
Ncellmasked = 0;          // -- number of masked cells = size(codim0masked,1)
codim0masked = [];
//
// "faces" are codim1 elements (interfaces between control volumes in FV method)
// - labelled as they appear during construction of cells:
//   left, right, down, top for each cells
//
// |--4--|--7--|--10-|--    3*Nx+1--|
// 1  1  2  2  5  3  8 =>    Nx  3*Nx-1
// |__3__|__6__|__9__|__      3*Nx__|
//
Nface = 0;          // -- number of faces = size(codim1,1)
codim1 = [];        // -- each row stores the (x,y) coordinates of the vertices
//                     i.e. v1_x,v1_y,v2_x,v2_y with [v1,v2] oriented from  
//                     bottom to top (vertical edges) or left to right (horizontal edges)
//
//codim0to1 = [];   // -- connectivity: non-zero codim0to1(i,j) stores i/j interface,
//            _antisymmetric_ matrix (when oriented), with morse reprsntatn A,B,face
codim0to1A = []; 
codim0to1B = [];
codim0to1NX = []; 
codim0to1NY = [];
codim0to1E = [];
//
Nghost = 0;
//
volume = (hx*hy);   // -- uniform here: all cells (elements) are similar rectangles

neighbours = 4;     // -- uniform here (except at boundaries: we artificially put -1)

diffusion = diffusion/volume;
source = source/volume; // CHECK

fboundary = list();
boundary = list();
periodic = list();
antiperiodic = list();
for i=1:4
  fboundary(i) = [];
  boundary(i) = [];
  periodic(i) = [];  
  antiperiodic(i) = [];
end

exec(mygeofile);

disp(Ncell,'Number of cells ');

disp(Nface,'Number of faces ');

// // // -- to build the skyline matrix with sparse Morse storage I,J,C=k
// // // codim0to1 = [codim0to1;zeros(Ncell,Ncell)];
// // // for k=1:Nface,
// // //   if((codim0to1I(k)>0)&(codim0to1J(k)>0)) // interior face
// // //     codim0to1(codim0to1A(k),codim0to1B(k)) = k;  // codim0to1C(k);
// // //     codim0to1(codim0to1B(k),codim0to1A(k)) = -k; //-codim0to1C(k);
// // //   end
// // // end

//******************************************************************************

// -- for error comp wrto analytical sol, only (fine) cell centers are needed // REMOVE
if(~onedy)
  lx = int(2^(levelref-mylevel)); // number of refinement in x
else
  lx = 1;
end
if(~onedx)
  ly = int(2^(levelref-mylevel)); // number of refinement in y
else
  ly = 1;
end
//
hxref = hx/lx; //hxref = Lx/(2^levelref);
hyref = hy/ly; //hyref = Ly/(2^levelref);
XXref = list();
YYref = list();
Npatterns = 0;
for llx=0:lx-1
  for lly=0:ly-1
    xx0 = (-(lx-1)/2+llx)*hxref; // -(hxref/2)*(lx-1)+llx*hxref;
    yy0 = (-(ly-1)/2+lly)*hyref; // -(hyref/2)*(ly-1)+lly*hyref;
    xref = [(x0+xx0):hx:(x0+xx0+(Nx-1)*hx)]; // -- periodic pattern !!
    yref = [(y0+yy0):hy:(y0+yy0+(Ny-1)*hy)]; // -- periodic pattern !!
    [YY,XX] = meshgrid(yref,xref);
    Npatterns = Npatterns+1;
    XXref(Npatterns) = matrix(XX,length(XX),1); // -- i=2;[XXref(i),YYref(i)]
    YYref(Npatterns) = matrix(YY,length(YY),1); // -- for check
  end
end
disp(Npatterns,'Number of patterns'); // lx*ly
if myhist & mylevel<levelmax then
  dlevel = levelmax-mylevel;
  mystart = 1;
  qq = [mystart:dlevel:mystart+2^levelmax];
  if ~onedx | ~onedy then
    while mystart<=(2^levelmax+1)*2^levelmax
      mystart = mystart + (2^levelmax+1)*dlevel;
      qq = [qq,[mystart:dlevel:mystart+2^levelmax]];
    end
  end
end // REMOVE

// ---------------------------------------------------------------------------
// ------------------------------ INITIALIZATION -----------------------------
// ---------------------------------------------------------------------------

clear qini q0 q1 qtp // Pdissipation

qini(1:Ncell,1:dimsys) = qana(0.,codim0(:,1),codim0(:,2)); // -- initial value for the unknown 

qini_todraw = list();
for i=1:dimsys
  qini_todraw(i) = matrix(state2plot(i,qini),[Nx,Ny]);
end

if sum(myplot(i))~=0 then clf();f0=gcf(); end

for i=1:dimsys
  if (myplot(i)) then
      subplot(1,sum(myplot*1.),sum(myplot(1:i)*1.));
    if (onedy) then
      ord = x0:hy:Ly;
      plot2d(ord,qini_todraw(i)(1,:),...
        sty(i),rect=[0,mymin(i),Ly,mymax(i)],leg=leg(i));
    elseif (onedx) then
      abc = y0:hx:Lx;
      plot2d(abc,qini_todraw(i)(:,1),...
        sty(i),rect=[0,mymin(i),Lx,mymax(i)],leg=leg(i));
    else
      abc = x0:hx:Lx;
      ord = y0:hy:Ly;
      // f=gcf();f.color_map=hotcolormap(10) // get handle of the parent figure
      xset("colormap",hotcolormap(10));
      // xset("colormap",jetcolormap(64))
      // Sgrayplot(abc,ord,qini_todraw(i),strf="041"); // smoothed	  
      plot3d1(abc,ord,qini_todraw(i)); // zoom_rect([0,0,1,1]);
      //colorbar(mymin(i),mymax(i));
      //contour2d(abc,ord,qini_todraw(i),...
      //  0:0.1:.6,nax=[3,(Nx-1)/4+1,3,(Ny-1)/4+1],rect=[0,0,1,1]);
    end
  end
end // end for i=1:dimsys

q0 = qini;               // -- current value
q1 = qini;               // -- one-step-older value

if myhist & mylevel==levelmax
  qhist = list();
  l0l0 = list();
  l0l1 = list();
  l0l2 = list();
  l2l1 = list();
  l2l2 = list();
  for i=1:dimsys
    qhist(i) = [];
    l0l0(i) = [];
    l0l1(i) = [];
    l0l2(i) = [];
    l2l1(i) = [];
    l2l2(i) = [];
  end
end

// preprocessing below not useful a priori, just in case
if (~onedx) then
  if (BC(1)=='periodic') then
    // if (BC(3)!='periodic') then error
    vctp = ( q0(periodic(1),1:dimsys)+q0(periodic(3),1:dimsys) )/2.;
    q0(periodic(1),1:dimsys) = vctp;
    q0(periodic(3),1:dimsys) = vctp;
//   elseif (BC(1)=='antiperiodic') then
//     // if (BC(3)!='antiperiodic') then error
//     vctp = ( q0(antiperiodic(1),1:dimsys)+q0(antiperiodic(3),1:dimsys) )/2.;
//     q0(antiperiodic(1),1:dimsys) = vctp;
//     q0(antiperiodic(3),1:dimsys) = vctp;
  end
end
if (~onedy) then
  if (BC(2)=='periodic') then
    // if (BC(4)!='periodic') then error
    vctp = ( q0(periodic(2),1:dimsys)+q0(periodic(4),1:dimsys) )/2.;
    q0(periodic(2),1:dimsys) = vctp;
    q0(periodic(4),1:dimsys) = vctp;
//   elseif (BC(2)=='antiperiodic') then
//     // if (BC(4)!='antiperiodic') then error
//     vctp = ( q0(antiperiodic(2),1:dimsys)+q0(antiperiodic(4),1:dimsys) )/2.;
//     q0(antiperiodic(2),1:dimsys) = vctp;
//     q0(antiperiodic(4),1:dimsys) = vctp;
  end
end

t = Tinitial;            // -- current time
nt = 0;                  // -- current iteration number
dnt = 0;                 // -- shift iteration number with respect to finest mesh
lastsnap = 0;            // -- current number of snapshots

coef = zeros(Ncell,1);          // -- stores maximal 1D speed wave   
flux = zeros(Ncell,dimsys+1);   // -- last col = "entropy" flux

// ---------------------------------------------------------------------------
// ------------------------------ COMPUTATIONS -------------------------------
// ---------------------------------------------------------------------------

while nt<Nsteps & t<Tfinal-%eps*10, // TIME LOOP //////////////////

  coef = coef*0;
  flux = flux*0;
  iface = 0;
  while iface<Nface, // FV : FACE LOOP FOR FLUXES COMPUTATION //////////////////
    iface = iface+1;
    
    // ** local change of frame (recall Galilean invariance)
    NX = codim0to1NX(iface);
    NY = codim0to1NY(iface);

    icell1 = codim0to1A(iface);   // shall receive a "left-type" flux
    icell2 = codim0to1B(iface);   // shall receive a"right-type" flux

    if (icell1<0)&(icell2>0) then // icell1 ghost
      // icell2 is a boundary cell, icell1 is fictitious    
      bndryflag = 1; 
      state2 = q0(icell2,:);      
      lambda1 = 1;
      lambda2 = 1;
      if (NY==0) then // vertical edge
        mytype = codim0(icell2,4);
        if (BC(mytype)=='periodic') then
          if (mytype==4)&(codim0(icell2+(Nx-1),4)~=2) then abort end
          state1 = q0(icell2+(Nx-1)-1,:); 
        elseif (BC(mytype)=='antiperiodic') then
          //if(icell2+Nx<=Ncell) then state1 = q0(icell2+Nx,:); else state1 = q0(icell2,:); end
          if(icell2-Nx>=1) then state1 = q0(icell2-Nx,:); 
            if (mytype==4)&(codim0(icell2-Nx,4)~=4) then abort end, 
            //state1 = q0(icell2,:); if state2~=state1 then abort end, 
          else 
            state1 = q0(icell2,:); 
          end
        else
          state1 = myBC(mytype,state2,codim0(icell2,1:2));
        end
      elseif (NX==0) then // horizontal edge 
        mytype = codim0(icell2,3);
        if (BC(mytype)=='periodic') then
          if (mytype==1)&(codim0(icell2+Nx*(Ny-1),3)~=3) then abort end
          state1 = q0(icell2+Nx*(Ny-2),:); 
        elseif (BC(mytype)=='antiperiodic') then
          //if(icell2+1<=Nx) then state1 = q0(icell2+1,:); else state1 = q0(icell2,:); end
          if(icell2-1>=1) then 
            state1 = q0(icell2-1,:); 
            if (mytype==1)&(codim0(icell2-1,3)~=1) then abort end, 
            //state1 = q0(icell2,:); if state2~=state1 then abort end, 
          else 
            state1 = q0(icell2,:); 
          end
        else
          state1 = myBC(mytype,state2,codim0(icell2,1:2));
        end
      else
        warning("error icell1");abort
      end	    
      if mymask(codim0(icell2,1),codim0(icell2,2)) then icell2=-Nghost-1; end 
    elseif (icell2<0)&(icell1>0) then // icell2 ghost
     // icell1 is a boundary cell, icell2 is fictitious
      bndryflag = 2;       
      state1 = q0(icell1,:); 
      lambda1 = 1;
      lambda2 = 1;
      if (NY==0) then // vertical edge
        mytype = codim0(icell1,4);
        if (BC(mytype)=='periodic') then
          if (mytype==2)&(codim0(icell1-(Nx-1),4)~=4) then abort end
          state2 = q0(icell1-(Nx-1)+1,:); 
        elseif (BC(mytype)=='antiperiodic') then
          //if(icell1-Nx>=1) then state2 = q0(icell1-Nx,:); else state2 = q0(icell1,:); end
          if(icell1+Nx<=Ncell) then state2 = q0(icell1+Nx,:);
            if (mytype==2)&(codim0(icell1+Nx,4)~=2) then abort end, 
            //state2 = q0(icell1,:); if state2~=state1 then abort end,
          else 
            state2 = q0(icell1,:); 
          end
        else
          state2 = myBC(mytype,state1,codim0(icell1,1:2));
        end
      elseif (NX==0) then // horizontal edge
        mytype = codim0(icell1,3);
        if (BC(mytype)=='periodic') then
          if (mytype==3)&(codim0(icell1-Nx*(Ny-1),3)~=1) then abort end
          state2 = q0(icell1-Nx*(Ny-2),:); 
        elseif (BC(mytype)=='antiperiodic') then
          //if(icell1-1>=Ncell-Nx) then state2 = q0(icell1-1,:); else state2 = q0(icell1,:); end
          if(icell1+1<=Ncell) then state2 = q0(icell1+1,:); 
            if (mytype==3)&(codim0(icell1+1,3)~=3) then abort end, 
            //state2 = q0(icell1,:); if state2~=state1 then abort end,
          else 
            state2 = q0(icell1,:); 
          end
        else
          state2 = myBC(mytype,state1,codim0(icell1,1:2));
        end
      else
        warning("error icell2");abort
      end
      if mymask(codim0(icell1,1),codim0(icell1,2)) then icell1=-Nghost-1; end 	          
    elseif (icell2>0)&(icell1>0)
    //   if (codim0(icell2,1)>codim0(icell1,1)|codim0(icell2,2)>codim0(icell1,2)) then    
      bndryflag = 0;
      // ** case without space-varying fluxes: lambda1 = 1;
      lambda1 = ( fux(t,codim0(icell1,1),codim0(icell1,2))*codim0to1NX(iface) + ...
		  fuy(t,codim0(icell1,1),codim0(icell1,2))*codim0to1NY(iface) );
      // ** case without space-varying fluxes: lambda2 = 1;   
      lambda2 = ( fux(t,codim0(icell2,1),codim0(icell2,2))*codim0to1NX(iface) + ...
		  fuy(t,codim0(icell2,1),codim0(icell2,2))*codim0to1NY(iface) );
      if ~mymask(codim0(icell2,1),codim0(icell2,2)) & mymask(codim0(icell1,1),codim0(icell1,2)) then  
        state2 = q0(icell2,:); 
        icell1 = -Nghost-1;
	state1 = myBC(0,state2,codim0(icell2,1:2));
      elseif mymask(codim0(icell2,1),codim0(icell2,2)) & ~mymask(codim0(icell1,1),codim0(icell1,2)) then
        state1 = q0(icell1,:); 
        icell2 = -Nghost-1;
	state2 = myBC(0,state1,codim0(icell1,1:2));
      elseif ~mymask(codim0(icell2,1),codim0(icell2,2)) & ~mymask(codim0(icell1,1),codim0(icell1,2)) then
        state1 = q0(icell1,:); 
        state2 = q0(icell2,:);
      else 
        icell1 = -Nghost-1;
        icell2 = -Nghost-2;
      end
    //   else warning("bad orientation assumed");
    //   end
    //else (icell2<0)&(icell1<0)
    end
    
    if (icell2>0)|(icell1>0) then 
      [leftf,rightf,lambda] = RIEMANN(lambda1,state1,lambda2,state2,NX,NY,bndryflag);
    end //else (icell2<0)&(icell1<0) // [h1,Q1,T1,X1,Y1,O1,Z1],  [h2,Q2,T2,X2,Y2,O2,Z2]
    
    if icell1>0 then // icell1 receives a "left-type" flux
      flux(icell1,1:dimsys+1) = flux(icell1,1:dimsys+1) + codim0to1E(iface)*leftf;
      if (CFLhard) then // -- "left"-type = -"right"-type in conservative case **
        coef(icell1) = max(coef(icell1),lambda);
      else // -- guermond: not stringent enough for entropy dissipation
        coef(icell1) = coef(icell1) + lambda*codim0to1E(iface)/neighbours;
      end
    end
    
    if icell2>0 then // icell2 receives a"right-type" flux
      flux(icell2,1:dimsys+1) = flux(icell2,1:dimsys+1) + codim0to1E(iface)*rightf;
      if (CFLhard) then // -- "left"-type = -"right"-type in conservative case **
	  coef(icell2) = max(coef(icell2),lambda);
      else // -- guermond: not stringent enough for entropy dissipation
	  coef(icell2) = coef(icell2) + lambda*codim0to1E(iface)/neighbours;
      end
    end

  end // FACE LOOP : FLUXES COMPUTATION ////////////////////// while iface<Nface  

  nt = nt+1;

  if (CFLhard) then // 1/2 CFL constraint !!
    k = min(Tfinal-t+%eps,volume/((hx*(~onedx)+hy*(~onedy))*2*max(coef)));
  else // *** guermond: not stringent enough for entropy dissip => 2 versions  
    k = min(Tfinal-t+%eps,volume/max(coef)); // 1/2 CFL also possible...
  end 
  
  if(tstpcst & levelmax>mylevel)
    if (k<TIME(nt,1)-t) then
      warning('Cannot use same time step as the finest grid');
    else
      k = TIME(nt,1)-t;
    end
    if dnt<>0 then 
      warning('error dnt');
    end    
  end
  
  t = t+k;

  q1 = q0;

  qtp = q1 + (k/volume)*flux(:,1:dimsys); // ** a piori faster than for loop in scilab
  
  // -- free-energy (entropy) dissipation after projection step (P) --  
  if globaltest & mylevel==levelmax then 
    F1 = myentropy2(q1(:,1:dimsys));
    F0 = myentropy2(qtp(:,1:dimsys));  
    Pdissipation(:,nt) = ( F0-F1-(k/volume)*flux(:,dimsys+1) );
    Pdissipationrel(:,nt) = Pdissipation(:,nt)./F1;
  end
  
  q0 = qtp;
  
// if viscosity>0 then
//   AA = sparse( speye(Nx*Ny,Nx*Ny)+viscosity*(k/volume)*diag(q0(:,1))*diffusion );
//   bb = sparse( qtp(:,2)+viscosity*(k/volume)*(q0(:,1).*(source*[0;0;1;0])) );
//   q0(:,2) = AA\bb;
// end

  if myhist & mylevel==levelmax then
    for i=1:dimsys
      qhist(i)(:,nt) = q0(:,i);
    end
  end
  
  ////////////////////////////////////////////////////////////////////////////////
  // -- POSTPROCESSING --  ///////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  
  TIME(nt+dnt,levelmax-mylevel+1) = t;
  
  NORM(nt+dnt,levelmax-mylevel+1) = sum(abs(q0-q1)); // ,'r'

  ddnt = 0;
  if( mylevel<levelmax )
    while t<Tfinal & TIME(nt+dnt+ddnt+1,1)<t
      ddnt = ddnt+1;
      TIME(nt+dnt+ddnt,levelmax-mylevel+1) = TIME(nt+dnt,levelmax-mylevel+1);
      NORM(nt+dnt+ddnt,levelmax-mylevel+1) = NORM(nt+dnt,levelmax-mylevel+1);
    end
    if(tstpcst & ddnt>1) 
      warning("error option time-step cst");
    end
  end  
  
  // -- OUTPUT: ERROR COMPUTATION, FIGURES
  
  qt = zeros(Nx*Ny,dimsys);
  for mycnt = 1:Npatterns;
    qt = qt + qana(t,XXref(mycnt),YYref(mycnt));
  end
  qt = qt/Npatterns; // -- weak (L^1, L^2 norm) type error: FV, FE
  
  nn = 1; // -- Npatterns = 2^(levelref+1) possibilities 
  qtt = qana(t,XXref(nn),YYref(nn)); // -- strong (pointwise) error: FD, or for plot
  
  q0_todraw = list();
  
  qref_todraw = list();

  if sum(myplot(i))~=0 then clf(0); f0=gcf(); end

  shoot = %F;
  if (t-lastsnap*period>=period) then // -- FOR PDE, GIFS ANIMATED
    lastsnap = lastsnap+1;
    shoot = %T;
    nlastsnap = string(lastsnap);
    if (lastsnap<10) then
      nlastsnap = "0"+nlastsnap;
    end
  end // -- end of if (t-lastsnap*period>=period) 
  
  if globaltest then
    disp(string(nt)+": "+string(t)+" (+"+string(k)+") dissip<="+string(max(Pdissipationrel(:,nt))));
  else
    disp(string(nt)+": "+string(t)+" (+"+string(k)+")");
  end
  
  disp("delta :"+string(NORM(nt+dnt,levelmax-mylevel+1)));
  
  for i=1:dimsys
    
    // -- to check conservativity : preservation of mass -- /////////////////////// 
    MASS(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel+1) = sum(q0(:,i)*volume);
      
    if(errorcomp(i)) // -- to check convergence num : rate of error decrease -- //////
      L0L0(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel+1) = volume*norm(q0(:,i)-qtt(:,i),%inf);
      L0L1(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel+1) = volume*norm(q0(:,i)-qt(:,i),1);
      L0L2(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel+1) = sqrt(volume*norm(q0(:,i)-qt(:,i),2));
      L2L1(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel+1) = (k/(ddnt+1))* ...
        (L0L1(i)(nt+dnt,levelmax-mylevel+1)^2);
      L2L2(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel+1) = (k/(ddnt+1))* ...
        L0L2(i)(nt+dnt,levelmax-mylevel+1)^2;
    end
    
    if myhist & mylevel<levelmax then // to check ORDER of scheme
      ql = qhist(i)(qq,nt+dnt);
      matrix(Ql,(2^mylevel)^2)/2^(levelmax-mylevel+1);
      l0l1(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel) = volume*norm(q0(:,i)-ql,1);
      l0l2(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel) = sqrt(volume*norm(q0(:,i)-ql,2));
      l2l1(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel) = (k/(ddnt+1))* ...
        l0l1(i)(nt+dnt,levelmax-mylevel)^2;
      l2l2(i)(nt+dnt:nt+dnt+ddnt,levelmax-mylevel) = (k/(ddnt+1))* ...
        l0l2(i)(nt+dnt,levelmax-mylevel)^2;
    end
   
    q0_todraw(i) = matrix(state2plot(i,q0),[Nx,Ny]);
    
    qref_todraw(i) = matrix(state2plot(i,qtt),[Nx,Ny]); // matrix(qt(:,i),[Nx,Ny]); 
    
    // -- to show monotonicity : preservation of invariant domains -- ////////////
    disp(leg(i)+" in ("+string(min(q0_todraw(i)))+";"+string(max(q0_todraw(i)))+")")

    if (myplot(i))&( (modulo(nt,plotevery)==0)|(nt==Nsteps)|(t>=Tfinal-%eps*10) ) then 

      currentfigurehandle = scf(0);
      
// // subplot(1,sum(myplot),sum(myplot(1:i))); // cfl(currentfigurehandle);
// // or
// if (sum(myplot)<=2) then
//   subplot(1,sum(myplot),sum(myplot(1:i))); // cfl(currentfigurehandle);
// elseif (i<4) then // add a flag for intensive var plot rather than extensive ??
//   subplot(2,sum(myplot)/2,sum(myplot(1:i))); // cfl(currentfigurehandle);
// else      
//   subplot(2,sum(myplot)/2,1.*(myplot(dimsys+1))+sum(myplot(1:i))); // cfl(currentfigurehandle);
// end  
      
      if (onedy) then
        plot2d([ord'],...
                [q0_todraw(i)(1,:)' qini_todraw(i)(1,:)'],...
                sty(i),leg=leg(i)+"@ini",rect=[0,mymin(i),Ly,mymax(i)]);
//                 [q0_todraw(i)(1,:)' qini_todraw(i)(1,:)' qref_todraw(i)(1,:)'],...
//                 sty(i),leg=leg(i)+"@ini@ref",rect=[0,mymin(i),Ly,mymax(i)]);
      //legends([leg(i);leg(i)+string("ini");leg(i)+string("ref")],sty(i),opt="lr")
        sleep(mysleep) // --  formerly xpause(100); in milliseconds      
      elseif (onedx) then
        plot2d([abc'],...
                [q0_todraw(i)(:,1) qini_todraw(i)(:,1)],...
                sty(i),leg=leg(i)+"@ini",rect=[0,mymin(i),Lx,mymax(i)]);
//                 [q0_todraw(i)(:,1) qini_todraw(i)(:,1) qref_todraw(i)(:,1)],...
//                 sty(i),leg=leg(i)+"@ini@ref",rect=[0,mymin(i),Lx,mymax(i)]);
        //zoom_rect([0,mymin(i),Lx,mymax(i)]);
        //legends([leg(i);leg(i)+string("ini");leg(i)+string("ref")],sty(i),opt="lr")
        //legends([string("256");string("128");string("64")],[1 2 3],opt="ul"); 
      else
        if mycontour then 
          xset("fpf"," ");
          // trick 1: this implies that the level numbers are not drawn on the level curves
          //ff.color_map=jetcolormap(7);
          xset("colormap",jetcolormap(8))
          // trick 2: to be able to put the legend on the right without interfering
          //          with the level curves use rect with a xmax value large enough
          //contour2d(abc,ord,q0_todraw(i),-0.2:0.2:1.2,nax=[7,(Nx-1)/8+1,7,(Ny-1)/8+1],rect=[0,0,1,1]);
          // trick 3: use legends (note that the more practical legend function
          //          will not work as soon as one of the level is formed by 2 curves)
          if withini then contour2d(abc,ord,qini_todraw(i),...
            0:0.1:.6,nax=[3,(Nx-1)/4+1,3,(Ny-1)/4+1],rect=[0,0,1,1]); end
          legends(string(0:0.1:.5),1:6,"lr"); 
          // legends(string(-0.2:0.2:1.2),1:8,"lr");
          xgrid(2);
          contour2d(abc,ord,q0_todraw(i),...
            0:0.1:.6,nax=[3,(Nx-1)/4+1,3,(Ny-1)/4+1],rect=[0,0,1,1]);
          //xs2png(0,myfilename+string(mylevel)+"nt"+nlastsnap+"contour"+leg(i)+".png") 
          //xs2jpg(0,myfilename+string(mylevel)+"nt"+nlastsnap+"contour"+leg(i)+".jpg") 
          //xs2eps(0,myfilename+string(mylevel)+"nt"+nlastsnap+"contour"+leg(i)+".eps") 
        else
          f=gcf(); //get the handle of the parent figure
          f.color_map=hotcolormap(10); // xset("colormap",hotcolormap(10))
          // xset("colormap",jetcolormap(64))
          // Sgrayplot(abc,ord,q0_todraw(i),strf="041"); // smoothed	  
          // colorbar(min(qini_todraw(i)),max(qini_todraw(i)))
          if withini then plot3d1(abc,ord,qini_todraw(i)); end
          //plot3d1(abc,ord,qref_todraw(i)); 
          plot3d1(abc,ord,q0_todraw(i)); // colorbar(min(q0_todraw(i)),max(q0_todraw(i)));
          // zoom_rect([0,0,1,1]);
          //xs2png(0,myfilename+string(mylevel)+"nt"+nlastsnap+"surface"+leg(i)+".png") 
          //xs2jpg(0,myfilename+string(mylevel)+"nt"+nlastsnap+"surface"+leg(i)+".jpg") 
          //xs2eps(0,myfilename+string(mylevel)+"nt"+nlastsnap+"surface"+leg(i)+".eps") 
        end
        xtitle(leg(i)); // show_window();
      end

    end // -- end of if (myplot(i))
    
  end // -- end of for i=1:dimsys
  
  if ( (modulo(nt,plotevery)==0)|(nt==Nsteps)|(t>=Tfinal-%eps*10) ) then 
  
    if myplot(dimsys+1) then
      subplot(2,sum(myplot)/2,4); // cfl(currentfigurehandle);
      if (onedy) then
	plot2d(ord',q0_todraw(4)(1,:).*q0_todraw(5)(1,:).*q0_todraw(7)(1,:),...
	      2,leg="det",rect=[0,0,Ly,10]);
	sleep(mysleep) // --  formerly xpause(100); in milliseconds      
      elseif (onedx) then
	plot2d(abc',q0_todraw(4)(:,1).*q0_todraw(5)(:,1).*q0_todraw(7)(:,1) ,...
	      2,leg="det",rect=[0,0,Lx,10]);
      else // plot the vector field
	//xset("colormap",hotcolormap(20))
	xset("colormap",jetcolormap(20))
	champ1(abc,ord,q0_todraw(2),q0_todraw(3));
//	champ1(abc,ord,q0_todraw(2),q0_todraw(3),rect=[0,0,1,1]);
//       colorbar(sqrt(min(q0_todraw(2).^2+q0_todraw(3).^2)),...
// 	       sqrt(max(q0_todraw(2).^2+q0_todraw(3).^2)))
	// xgrid(1)//grille
	//       cmap=jetcolormap(64)
	//       f=gcf()//choix table couleurs
	//       f.color_map=cmap
      end
    end

    fprintfMat("./"+myfilename+"time.dat",TIME);
    fprintfMat("./"+myfilename+"norm.dat",NORM);
    fprintfMat("./"+myfilename+"q0.txt",q0);
    fprintfMat("./"+myfilename+"q1.txt",q1);
    fprintfMat("./"+myfilename+"X.txt",codim0(:,1:2));

  end
    
  if sum(myplot)~=0  then  
    if (onedx|onedy) then
      if shoot & mygif then 
        xs2gif(0,myfilename+string(mylevel)+"nt"+nlastsnap+"1D.gif") 
      end
      if shoot & mypdf then 
        xs2pdf(0,myfilename+string(mylevel)+"nt"+nlastsnap+"1D.pdf") 
      end
    else
      if mycontour then
        if shoot & mygif then 
          xs2gif(0,myfilename+string(mylevel)+"nt"+nlastsnap+"contour.gif") 
        end
        if shoot & mypdf then 
          xs2pdf(0,myfilename+string(mylevel)+"nt"+nlastsnap+"contour.pdf") 
        end
      else
        if shoot & mygif then 
          xs2gif(0,myfilename+string(mylevel)+"nt"+nlastsnap+"surface.gif") 
        end
        if shoot & mypdf then 
          //xs2pdf(0,myfilename+string(mylevel)+"nt"+nlastsnap+"surface.pdf") 
          xs2png(0,myfilename+string(mylevel)+"nt"+nlastsnap+"surface.png") 
        end
      end
    end
    sleep(mysleep) // --  formerly xpause(100); in milliseconds      
  end  
  
  dnt = dnt+ddnt; // update nb iter shift wrto levelmax
  
end //  TIME LOOP ////////////////////////////// while t<Tfinal-%eps & nt<Nsteps

ddnt = 0;
while(nt+dnt+ddnt<size(TIME,1)) // TO POLISH POST-PROCESSING: NOT USELESS !!
  ddnt = ddnt+1;
  TIME(nt+dnt+ddnt,levelmax-mylevel+1) = TIME(nt+dnt,levelmax-mylevel+1); // t ??
  NORM(nt+dnt+ddnt,levelmax-mylevel+1) = NORM(nt+dnt,levelmax-mylevel+1);
  for i=1:dimsys
    MASS(i)(nt+dnt+ddnt,levelmax-mylevel+1) = MASS(i)(nt+dnt,levelmax-mylevel+1);
    L0L0(i)(nt+dnt+ddnt,levelmax-mylevel+1) = L0L0(i)(nt+dnt,levelmax-mylevel+1);
    L0L1(i)(nt+dnt+ddnt,levelmax-mylevel+1) = L0L1(i)(nt+dnt,levelmax-mylevel+1);
    L0L2(i)(nt+dnt+ddnt,levelmax-mylevel+1) = L0L2(i)(nt+dnt,levelmax-mylevel+1);
    L2L1(i)(nt+dnt+ddnt,levelmax-mylevel+1) = L2L1(i)(nt+dnt,levelmax-mylevel+1);
    L2L2(i)(nt+dnt+ddnt,levelmax-mylevel+1) = L2L2(i)(nt+dnt,levelmax-mylevel+1);
  end
end

if mylevel==levelmax & t<Tfinal then
  if nt==Nsteps then
    warning("Final time changed due to Nsteps on finest grid")
    Tfinal = t;
  else
    error("Problem final time on finest grid")
  end
end

end ////////////////////////////////////////// for mylevel=levelmax:-1:levelmin

// -----------------------------------------------------------------------------
// --------------------------------- OUTPUTS ----------------------------------
// -----------------------------------------------------------------------------

for i=1:dimsys

  if (massshow(i)) then
  
    ff = scf(i);

    plot2d(TIME,MASS(i)); title("Mass "+leg(i))

    fprintfMat("./MASS"+leg(i)+".dat",MASS(i));
    
  end

  if(errorcomp(i)) then
    
    fprintfMat("./L0L1"+leg(i)+".dat",L0L1(i)); // -- time-inst L1 conv
    fprintfMat("./L0L2"+leg(i)+".dat",L0L2(i));
    fprintfMat("./L2L1"+leg(i)+".dat",L2L1(i)); // -- l2-time L1 conv
    fprintfMat("./L2L2"+leg(i)+".dat",L2L2(i)); // -- l2-time L2 conv

    for j=1:size(L0L1(i),1) // -- to compute speed of convergence in various norms
        [a,b,sig]=reglin([levelmax:-1:levelmin]*log(2),log(L0L1(i)(j,:)));
        aL0L1(i,j)=a; bL0L1(i,j)=b; sL0L1(i,j)=sig;
        [a,b,sig]=reglin([levelmax:-1:levelmin]*log(2),log(L0L2(i)(j,:)));
        aL0L2(i,j)=a; bL0L2(i,j)=b; sL0L2(i,j)=sig;
        //   [a,b,sig]=reglin([levelmax:-1:levelmin]*log(2),log(L2L1(i)(j,1:$-1)));
        //   aL2L1(i,j)=a; bL2L1(i,j)=b; sL2L1(i,j)=sig;
        //   [a,b,sig]=reglin([levelmax:-1:levelmin]*log(2),log(L2L2(i)(j,1:$-1)));
        //   aL2L2(i,j)=a; bL2L2(i,j)=b; sL2L2(i,j)=sig;
    end  // -- in a restricted time domain, space ?

    ff = scf(10+i); 
    
    clf();
    plot2d(TIME,L0L1(i)); // title(L0L1)
    xs2pdf(gcf(),'L0L1'+leg(i)+'.pdf'); // xs2png(i,'L0L1.png');
    clf(); 
    plot2d(TIME(:,1),[aL0L1(i,:)',sL0L1(i,:)']); 
    xs2pdf(gcf(),'L0L1confidencespeed'+leg(i)+'.pdf');
    clf(); 
    plot2d(TIME(:,1),aL0L1(i,:));
    xs2pdf(gcf(),'L0L1speed'+leg(i)+'.pdf');

    ff = scf(20+i);
    
    clf(); 
    plot2d(TIME,L0L2(i)); // zoom_rect([%eps,%eps,max(TIME(i)),max(L0L2(i))]);
    xs2pdf(gcf(),'L0L2'+leg(i)+'.pdf'); //xs2png(i,'L0L2.png');
    clf(); 
    plot2d(TIME(:,1),[aL0L2(i,:)',sL0L2(i,:)']); 
    xs2pdf(gcf(),'L0L2confidencespeed'+leg(i)+'.pdf');
    clf(); 
    plot2d(TIME(:,1),aL0L2(i,:)');
    xs2pdf(gcf(),'L0L2speed'+leg(i)+'.pdf');

    // clf(); plot2d(TIME(1:$-1,:),sqrt(cumsum(L2L1(i)(1:$-1,:),'c')));xs2pdf(1,'L2L1.pdf');
    // xs2pdf(gcf(),'L2L1.pdf'); //xs2png(i,'L2L1.png');
    // clf(); plot2d(TIME(:,1),[aL2L1(i,:)',sL2L1(i,:)']); 
    // clf(); plot2d(TIME(:,1),aL2L1(i,:)');
    // xs2pdf(gcf(),'L2L1speed.pdf');
    // 
    // clf(); plot2d(TIME(1:$-1,:),sqrt(cumsum(L2L2(i)(1:$-1,:),'c')));xs2pdf(1,'L2L2.pdf');
    // xs2pdf(i,'L2L2.pdf'); //xs2png(i,'L2L2.png');
    // clf(); plot2d(TIME(:,1),[aL2L2(i,:)',sL2L2(i,:)']); 
    // clf(); plot2d(TIME(:,1),aL2L2(i,:)');
    // xs2pdf(gcf(),'L2L2speed.pdf');

    for j=1:size(L0L0(i),1) // -- to compute speed of convergence in various norms
        [a,b,sig]=reglin([levelmax:-1:levelmin]*log(2),log(L0L0(i)(j,:)));
        aL0L0(i,j)=a; bL0L0(i,j)=b; sL0L0(i,j)=sig;
    end

    ff = scf(30+i);

    clf(); 
    plot2d(TIME,L0L0(i));
    xs2pdf(gcf(),'L0L0'+leg(i)+'.pdf');
    clf(); 
    plot2d(TIME(:,1),[aL0L0(i,:)',sL0L0(i,:)']); 
    xs2pdf(gcf(),'L0L0confidencespeed'+leg(i)+'.pdf');
    clf(); 
    plot2d(TIME(:,1),aL0L0(i,:)');
    xs2pdf(gcf(),'L0L0speed'+leg(i)+'.pdf');

    if myhist then
        for j=1:size(l0l1,1) // -- to compute order of the method in various norms: always better than speed !!
            [a,b,sig]=reglin([levelmax-1:-1:levelmin]*log(2),log(l0l1(j,:)));
            al0l1(i,j)=a; bl0l1(i,j)=b; sl0l1(i,j)=sig;
            [a,b,sig]=reglin([levelmax-1:-1:levelmin]*log(2),log(l0l2(j,:)));
            al0l2(i,j)=a; bl0l2(i,j)=b; sl0l2(i,j)=sig;
            //   [a,b,sig]=reglin([levelmax-1:-1:levelmin]*log(2),log(L2L1(j,:)));
            //   al2l1(i,j)=a; bl2l1(i,j)=b; sl2l1(i,j)=sig;
            //   [a,b,sig]=reglin([levelmax-1:-1:levelmin]*log(2),log(L2L2(j,:)));
            //   al2l2(i,j)=a; bl2l2(i,j)=b; sl2l2(i,j)=sig;
        end
    end

  end

end
