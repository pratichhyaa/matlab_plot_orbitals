clearvars
clc
clear global;
close all force;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  March 10 2026 
%%%%%%%%%%
%%%%%%%%%%  The program should read a Dyson or parent orbital through 
%%%%%%%%%%  qcout_Dyson_1.m
%%%%%%%%%%  Use Ylm_Reg for the pure harmonics (normalization is now correct)
%%%%%%%%%%  Currrently the order is s, px py pz, d0,d1,d-1,d2,d-2 etc
%%%%%%%%%%  (Can replace with module qcout to read a QChem result)
%%%%%%%%%%
%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%  module  molden.m to read the input                   %%%%%%%%%%
%%%%%%%%%%  atms = #of atoms                                     %%%%%%%%%%
%%%%%%%%%%  A = array of atom positions (geometry)               %%%%%%%%%%
%%%%%%%%%%  SPDFG - numbers of each type of orbital per atom     %%%%%%%%%%
%%%%%%%%%%  CZ - # of primitives and angular momentum for each AO%%%%%%%%%%
%%%%%%%%%%  CCEX - exponents and contraciton coefficents for     %%%%%%%%%%
%%%%%%%%%%         Gaussian Primitives                           %%%%%%%%%%
%%%%%%%%%%  DL, DR - AO coefficients for left and right DOs      %%%%%%%%%%
%%%%%%%%%%  llim - highest L in basis set                        %%%%%%%%%%
%%%%%%%%%%  AOS - number of atomic orbitals in basis set         %%%%%%%%%%
%%%%%%%%%%  MO - orbital index (start at 1)                      %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MO=1;
[filename1,Atoms,atms,A,SPDFG,CZ,CCEX,DL,DR,llim,AOS]=molden(MO);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%  Parameters for the Orbital(s)                        %%%%%%%%%%
%%%%%%%%%%  ps = #o points on grid                               %%%%%%%%%%
%%%%%%%%%%  iso - isovalue to display the orbital                %%%%%%%%%%
%%%%%%%%%%  gridmax = max distance from origin                   %%%%%%%%%%
%%%%%%%%%%  Dystol: tolerance for basis function ceofficients    %%%%%%%%%%
%%%%%%%%%%  npts - add 1 tp ps to include center                 %%%%%%%%%%
%%%%%%%%%%  stepsz: gridpoint spacings                           %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps=50;
iso=0.001;
gridmax=20;
Dystol=1E-15;
tol=8;
npts=ps+1;                            %number of gridpoints (per dimension)
stepsz=2*gridmax/ps;                              %stepsize for integration
cen=ps/2+1;                                                    %grid center
x1=-gridmax:stepsz:gridmax;                      %3D (cubic) grid dimension
[X,Y,Z]=meshgrid(x1);
r(npts,npts,npts)=zeros;
r(:,:,:)=sqrt(X.^2+Y.^2+Z.^2);                      %(r values on the grid)
a0=0.5291172109;                                               %Bohr Radius
A=A/a0;            %convert atom coordinates from Angstroms to atomic units
% num_lm = (llim + 1)^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%          Load Regular Solid Harmonics for each atom             %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars S rmat
% S = zeros(npts,npts,npts,ps,num_lm,'single');
for a=1:atms
rmatsq(:,:,:,a)=(X-A(a,1)).^2+(Y-A(a,2)).^2+(Z-A(a,3)).^2;
rmat(:,:,:,a)=sqrt(rmatsq(:,:,:,a));
S(:,:,:,:,:,a)=Ylm_Reg(ps,llim,X-A(a,1),Y-A(a,2),Z-A(a,3),rmat(:,:,:,a));
% S = S + single(Ylm_Reg(ps, llim, X - A(a,1), Y - A(a,2), Z - A(a,3), rmat));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%         Translate the Input into Orbitals on grid               %%%%%
%%%%%                                                                 %%%%%
%%%%%  WF - (:,:,:,1) parent (Koopman's picture) or LDyson            %%%%%
%%%%%  WF - (:,:,:,2) parent RDyson            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars WF
  WF(npts,npts,npts,2)=zeros;
            
AO=0;
l=0;
count=0;
for a=1:atms
    j=0; %index for Gaussian primitives 
    k=0; %index for Dyson coefficients
    for row=1:size(CZ,1)
        if CZ(row,a,1)==0
            break
        end
        l=CZ(row,a,2);
        clearvars AOmat %prepare the radial part for the given AO of atom a
        AOmat(npts,npts,npts)=zeros;
        AOmatang(npts,npts,npts)=zeros;
        for i=1:CZ(row,a,1)
            j=j+1;
             AOmat(:,:,:)=AOmat(:,:,:)+CCEX(j,a,2)*sqrt((sqrt((CCEX(j,a,1)^3)/pi))*(2^(7/2))).*exp(-CCEX(j,a,1)*rmatsq(:,:,:,a));
        end
        for m=-l:l  %account for/include the different m components of each AO
            AO=AO+1;
               WF(:,:,:,1)=WF(:,:,:,1)+DL(AO,1)*AOmat.*S(:,:,:,l+1,m+llim+1,a);
               WF(:,:,:,2)=WF(:,:,:,2)+DR(AO,1)*AOmat.*S(:,:,:,l+1,m+llim+1,a);
clearvars AOmatang
            end
clearvars AOmat
    end

end

cm=cmaps(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% rotate the orbital so that x or y or z MF lies along Z LF    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0;
beta=pi/2;
gamma = pi/2;
 Rmat=[(cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)),...
         (-cos(gamma)*sin(alpha)-cos(alpha)*cos(beta)*sin(gamma)),cos(alpha)*sin(beta);...
         (cos(alpha)*sin(gamma)+cos(beta)*cos(gamma)*sin(alpha)),...
         (cos(alpha)*cos(gamma)-cos(beta)*sin(alpha)*sin(gamma)),sin(alpha)*sin(beta);...
         -cos(gamma)*sin(beta),sin(beta)*sin(gamma),cos(beta)];
   

xLF=X*Rmat(1,1)+Y*Rmat(1,2)+Z*Rmat(1,3);
yLF=X*Rmat(2,1)+Y*Rmat(2,2)+Z*Rmat(2,3);
zLF=X*Rmat(3,1)+Y*Rmat(3,2)+Z*Rmat(3,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%% Setting up orientation
%%%%%%%%%% (this is the time consuming bit. By defining the original
%%%%%%%%%% wavefunction within a spherical volume of the grid, linear
%%%%%%%%%% or cubic interpolation is possible, quicker than splines)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars WFr
WFr=zeros(npts,npts,npts,2);
WFr(:,:,:,1)=interp3(X,Y,Z,WF(:,:,:,1),xLF,yLF,zLF,'linear',0);
WFr(:,:,:,2)=interp3(X,Y,Z,WF(:,:,:,2),xLF,yLF,zLF,'linear',0);
 vecz(:,1)=transpose(Rmat)*[0;0;gridmax-1];
 vecy(:,1)=transpose(Rmat)*[0;gridmax-1;0];
 vecx(:,1)=transpose(Rmat)*[gridmax-1;0;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%% Display the results through module quickvis
%%%%%%%%%% quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,gridmax,me_et,LR,hal,orb)
%%%%%%%%%% options are me_et - methyl (1) or ethyl (2) (basically tells how
%%%%%%%%%%                     many atoms to display
%%%%%%%%%% LR - which Dyson orbital (1) = L, (2) = R, use 1 for HF orbital
%%%%%%%%%% hal - which alkyl halide (1) = Cl, (2) = Br, (3) = I
%%%%%%%%%% orb - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if atms==6
    me_et = 1;
    if Atoms{5,1}=="Cl"
        hal=1;
    end
    if Atoms{5,1}=="Br"
        hal=2;
    end
    if Atoms{5,1}=="I"
        hal=3;
    end

    
end



if atms==9
    me_et=2
    if Atoms{8,1}=="Cl"
        hal=1;
    end
    if Atoms{8,1}=="Br"
        hal=2;
    end
    if Atoms{8,1}=="I"
        hal=3;
    end   
end
% orb=1;
LR=1;
%  quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,25,me_et,1,hal,0,MO)
%  pause(5)
% quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,25,me_et,1,hal,1,MO)
% pause(10)
% quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,25,me_et,1,hal,2,MO)
 quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,50,me_et,LR,hal,1,MO);
% quickvis(filename1,X,Y,Z,WF,npts,cm,iso,A,alpha,beta,gamma,25,me_et,1,hal,2,MO)