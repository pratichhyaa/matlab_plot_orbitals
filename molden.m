function [filename1,Atoms,atms,Coords,SPDFG,CZ,CCEX,DL,DR,llim,AOS]=molden(MO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%  Program to Read in Dyson Orbitals from QChem .out    %%%%%%%%%%
%%%%%%%%%%  .out file needs to have the basis set printed        %%%%%%%%%%
%%%%%%%%%%  .out file needs to have the DO coefficeints for th   %%%%%%%%%%
%%%%%%%%%%       AO basis                                        %%%%%%%%%%
%%%%%%%%%%  .out needs the geometry specified in the $molecule   %%%%%%%%%%
%%%%%%%%%%       section (i.e. the calculation was SP            %%%%%%%%%%
%%%%%%%%%%  [Can handle user defined or std basis set input]     %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%        WORKING VERSION MARCH 27 2024                  %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clearvars -except MO
clc
filepath1 = 'C:\Users\rmlab.CHEMISTRY\Desktop\plotting_molden_file\'; 
filename1 = 'ccsd_eomee_ccsd_rhfref_singlets_1_1_nto.out';

% 2. Combine them into an absolute path to guarantee MATLAB finds it
fullfilepath = fullfile(filepath1, filename1);

fid = fopen(fullfilepath, 'rt');
% [filename1,filepath1]=uigetfile({'*.out','Text Files'},...
%   'Select Data File 1');
%    cd(filepath1);
% fid = fopen(filename1, 'rt');
   if fid==-1
    disp('Error')
   end
   clearvars CCEX      
disp('Reading Output File');
i=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%  Read the QChem output file or Molden file            %%%%%%%%%%
%%%%%%%%%% line by line then close file                          %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bs=1;
    while feof(fid)==0
        A{i,1}=fgetl(fid); 
        i=i+1;
    end
    fclose(fid);
    length=size(A,1);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%         Find section with atom coordinates            %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        i=0;
    while i<length
        i=i+1;
        Astr=regexp(string(A(i)),'\S+','match');
        if numel(Astr)==2
        compare=strcmp(Astr,["[Atoms]" "(Angs)"]);
            if compare==1
                 break
            end
        end
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%                Count number of atoms                  %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        i=i+1;
        atms=0;
    while i<length
        i=i+1;
        Astr=regexp(string(A(i)),'\S+','match');
        compare=strcmp(Astr,["[GTO]"]);
                if compare == 1
                break 
            end
        atms=atms+1;
    end
    atms=atms+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%            Read coordinates of each atom              %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Atoms=cell(atms,1); %array for atom names - may be redundant
    Coords=zeros(atms,3); %array for coordinates of each atom (row x,y,z)
    start=i-atms-1; %rewind to start of geometry section
    for j=1:atms
        k=start+j;
        Astr=regexp(string(A(k)),'\S+','match');
        Atoms{j,1}=Astr(1);
        Coords(j,1)=str2double((Astr(4)));
        Coords(j,2)=str2double((Astr(5)));
        Coords(j,3)=str2double((Astr(6)));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%     Set up count of AOs of each L for each atom       %%%%%%%%%%
%%%%%%%%%%     Currently up to G, could go higher by increasing  %%%%%%%%%%
%%%%%%%%%%     number in SPDFG, and letters in angs              %%%%%%%%%%
%%%%%%%%%%     Read Basis Set                                    %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SPDFG=zeros(atms,5);
    angs=["S","P","D","F","G"];
    at=1;
    b=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%               Find start of basis set                   %%%%%
        %%%%%     (Currently not set up to deal with an ecp)          %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i=i+2;
    l=1;
    m=1;
  %  i=i+2;
    llim=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%                   read basis set                        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while i<length
        Astr=regexp(string(A(i)),'\S+','match');
        if isempty(Astr)
            Astr=regexp(string(A(i+1)),'\S+','match');
            if(strcmp(Astr(1),"[MO]")==1)
            break
            end
            i=i+2;
            at=at+1;
            l=1;
            m=1;

                 Astr=regexp(string(A(i)),'\S+','match');   
        end
         if at==atms+1
                 break
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%    is this a header line if so, how many primitives?    %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if ismember((Astr(1)),angs)==1
            CZ(m,at,1)=str2double((Astr(2)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%    what is the angular momentum of this orbital?        %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        if(strcmp(Astr(1),"S")==1)
            SPDFG(at,1)=SPDFG(at,1)+1;
            CZ(m,at,2)=0;
            i=i+1;
            m=m+1;
        end    

        if(strcmp(Astr(1),"P")==1)
            SPDFG(at,2)=SPDFG(at,2)+1;
            CZ(m,at,2)=1;
            i=i+1;
            m=m+1;
            if llim<1 
                llim=1;
            end
        end    
        if(strcmp(Astr(1),"D")==1)
            SPDFG(at,3)=SPDFG(at,3)+1;
            i=i+1;
            CZ(m,at,2)=2;
            m=m+1;
            if llim<2 
                llim=2;
            end
        end    
        if(strcmp(Astr(1),"F")==1)
            SPDFG(at,4)=SPDFG(at,4)+1;
            i=i+1;
            CZ(m,at,2)=3;
            m=m+1;
            if llim<3 
                llim=3;
            end
        end
        if(strcmp(Astr(1),"G")==1)
            SPDFG(at,5)=SPDFG(at,5)+1;
            i=i+1;
            CZ(m,at,2)=4;
            m=m+1;
        if llim<4 
            llim=4;
        end
        end    
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  if not header row or **** then the following will be   %%%%%
        %%%%%  exponents followed by contraction coefficients         %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if sum(ismember((Astr(1)),angs))==0
        CCEX(l,at,1)=str2double((Astr(1))); 
        CCEX(l,at,2)=str2double((Astr(2))); 
        l=l+1;
        i=i+1;
        end
           if (strcmp(Astr(1),"Ene")==1)
           break     
        end
    end
    display('done')
   CZ;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%     Find and read the Dyson Orbitals                  %%%%%%%%%%
%%%%%%%%%%     Left first then right, in AO basis                %%%%%%%%%%
%%%%%%%%%%     Key word in file is 'Decomposition'               %%%%%%%%%%
%%%%%%%%%%     Currently only the first transition               %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 
compare=0;
    i=i+1;
    i=i+4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%             Calculate # of AOs in AO basis              %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    AOS=sum(sum(SPDFG(:,1)))+3*sum(sum(SPDFG(:,2)))+...
        5*sum(sum(SPDFG(:,3)))+7*sum(sum(SPDFG(:,4)))+...
        9*sum(sum(SPDFG(:,5)))
    i=i+(MO-1)*AOS;
 i=i+(MO-1)*4;
 for k=1:atms
 sz(k)=sum(sum(SPDFG(k,1)))+3*sum(sum(SPDFG(k,2)))+...
        5*sum(sum(SPDFG(k,3)))+7*sum(sum(SPDFG(k,4)))+...
        9*sum(sum(SPDFG(k,5)));
 end
 h=max(sz);
 DLarray=zeros(k,atms);
    for k=1:atms
       
    for j=1:sz(k)
        i=i+1;
        Astr=regexp(string(A(i)),'\S+','match');
%        if isempty(Astr)==0
            DLarray(j,k)=str2double((Astr(2)));
%        end
    end
    end
DLarray;

i=i+4;
    DR=zeros(AOS,1);
%     for j=1:AOS
%         i=i+1;
%         Astr=regexp(string(A(i)),'\S+','match');
% %        if isempty(Astr)==0
% 
%             DR(j)=str2double((Astr(1)));
% %        end
%     end
%     disp('Input File Read');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%            Reorder atoms and coordinates              %%%%%%%%%%
%%%%%%%%%%            (H, C, hal, I)                             %%%%%%%%%%
%%%%%%%%%%                                                       %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
disp('Reordering atoms for display');

reorr=zeros(1,atms);
ids={"H" "C" "Cl" "Br" "I"};
reorr(1,:)=1:atms;
for i=1:5
    for j=1:atms
compare=strcmp(Atoms{j,1},ids{1,i});
if compare==1
    reorr(1,j)=i;
end
    end
end
[B,I]=sort(reorr,2);

Atomsnew(:,:)=Atoms(I,:);
coordnew(:,:)=Coords(I,:);

SPDFGnew(:,:)=SPDFG(I,:);
CZnew(:,:,:)=CZ(:,I,:);
CCEXnew(:,:,:)=CCEX(:,I,:);
Coords=coordnew;
SPDFG=SPDFGnew;
CZ=CZnew;
CCEX=CCEXnew;
Atoms=Atomsnew;
DL=zeros(AOS,1);
i=1;
for k = 1:atms
for j = 1:sz(I(k))
DL(i,1)=DLarray(j,I(k));
i=i+1;
end
end

DL;   
clearvars coordnew SPDFGnew CZnew CCEXnew Atomsnew

end