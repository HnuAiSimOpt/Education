%==========================================================================%
%                  HEAD                                                    %
%              TestPlatform_Implicit                                       %
%          CaoYeHong                                                       %
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
    nes = length(nodes(:,1));               % ϵͳ��Ԫ��
    nns = length(gcoords(:,1));             % ϵͳ�ڵ���
    nds = nns * ndpn;                       % ϵͳ���ɶ���
   ndpe = ndpn * nnpe;                      % ��Ԫ���ɶ���
     ff = sparse(nds,1);
     kk = sparse(nds,nds);
   Disp = sparse(nds,1);
  index = sparse(ndpe,1);
%==========================================================================

% CONSTITUTIVE MODELS

switch Element_Type
    case 'CPS6'
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
