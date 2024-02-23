%==========================================================================%
%                  HEAD                                                    %
%              TestPlatform_Implicit                                       %
%          CAOYEHONG                                                       %
%              HUNAN UNIVERSITY.                                           %
%                                                                          %
%==========================================================================%

% START
clear;
close all;
clc;
format long
%==========================================================================

% LOAD DATA
addpath('Results_Model');
selection=selectElement();
[gcoords,nodes,nnpe,ndpn,ElasticM,PoissonR,BC,load,dimension,Element_Type,SecType,Thickness,Density] = loadinp(selection);
fprintf('Load Data Finished\n');
%==========================================================================

% PAPERING WORK
    nes = length(nodes(:,1));               % 系统单元数
    nns = length(gcoords(:,1));             % 系统节点数
    nds = nns * ndpn;                       % 系统自由度数
   ndpe = ndpn * nnpe;                      % 单元自由度数
     ff = sparse(nds,1);
     kk = sparse(nds,nds);
   Disp = sparse(nds,1);
  index = sparse(ndpe,1);
%==========================================================================

% CONSTITUTIVE MODELS

switch Element_Type
    case 'S4R'
        Cm = Thickness*fematiso(1,ElasticM,PoissonR);
        Cb = Cm*Thickness*Thickness/12;
    shearm = 0.5*ElasticM/(1.0+PoissonR);
      scof = 5/6;
        Cs = shearm*scof*Thickness*eye(2);
    case 'C3D4'
        Cm = fematiso(4,ElasticM,PoissonR);
    case 'CPS6'
        Cm = fematiso(1,ElasticM,PoissonR);
    case 'CPS4'
        Cm = fematiso(1,ElasticM,PoissonR);
    otherwise
        disp('NOT SUPPORTED YET\n');
end
%==========================================================================


% COMPUTE THE STIFFNESS MATRIX OF SYSTEM.
startT = cputime;
switch Element_Type
    case 'CPS6'
        kk=getK_CPS6(gcoords,nodes,nes,Cm,nnpe,ndpn,kk);
    case 'S4R'
        kk = getK_S4R(gcoords,nodes,nes,Cm,Cb,Cs,Thickness,nnpe,ndpn,kk);
    case 'CPS4'
        kk=getK_CPS4(gcoords,nodes,nes,Cm,nnpe,ndpn,kk);
    case 'C3D4'
        kk=getK_C3D4(gcoords,nodes,nes,Cm,nnpe,ndpn,kk);
    otherwise
        fprintf('This element:%s is not supported yet!\n',Element_Type);
end
%==========================================================================


% APPLY LOAD
for count = 1:length(load(:,1))
    ff(load(count,1))=load(count,2);
end
%==========================================================================


% APPLY BOUNDARY CONDITION
kk_f = kk;
[kk,ff] = loadbc(kk,BC,ff);
%==========================================================================

% SOLVE EQUATION
fprintf('LU decomposition......\n');
[LL,UU] = lu(kk);
fprintf('Solving......\n');
Utemp = LL\ff;      
Disp = UU\Utemp;
clear LL UU Utemp;
TIME = cputime-startT;
fprintf('Computation Finished\nTime Cost:%fs\n',TIME);
%==========================================================================


% OUTPUT DATA
femout(Disp,gcoords,nodes,nns,ndpn,Element_Type);
%=========================================================================%
%                                END CODE                                 %
%=========================================================================%
