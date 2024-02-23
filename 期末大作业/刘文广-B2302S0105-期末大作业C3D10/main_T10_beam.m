% C3D10
% Model File：3Dbeam.inp
%------------------------------------
%Variable descriptions
% k = element stiffness matrix
% kk = system stiffness matrix
% ff = system force vector
% disp = system nodal displacement vector
% gcoord = coordinate values of each node
% nodes = nodal connectivity of each element
% index = a vector containing systems dofs associated with element
% bcdof = a vector containing dofs associated with boundary conditions
% bcval = a vector containing boundary condition values associated with the dofs in bcdof
% ---------------------------------------------------------------------------------------
clear
format long
first_time=cputime;
%----------------------------------
% input data for control parameters
%----------------------------------
L=10;                                    % length of the problem in x-axis 
B=1;                                     % length of the problem in y-axis 
H=1;                                     % length of the problem in z-axis 
emodule=2.1e9;                           % elastic modulus
poisson=0.3;                             % Poisson's ratio
p0=-1;                                % magnitude of the external load

nnel=10;                                 % number of node per element
ndof=3;                                  % number of dofs per node 
edof=nnel*ndof;                          % number of dofs per element 
%--------------------------
% read the associated datas
%--------------------------
fid_inp=fopen('3Dbeam.inp','r');         % open the input file
linestr=write_echo_file(fid_inp);
nnode=str2num(linestr);                  % total number of node 
linestr=write_echo_file(fid_inp);
nel=str2num(linestr);                    % total number of element
sdof=nnode*ndof;                         % total number of dofs 
%-----------------------------
% read nodal coordiante values
%-----------------------------
gcoord=zeros(nnode,3);
for inode=1:nnode
    clear tmp
    linestr=write_echo_file(fid_inp);
    tmp=str2num(linestr);
    gcoord(inode,1)=tmp(2);
    gcoord(inode,2)=tmp(3);
    gcoord(inode,3)=tmp(4);
end
%-----------------------------------------
% read nodal connectivity for each element
%-----------------------------------------
nodes=zeros(nel,nnel);
for iel=1:nel
    clear tmp
    clear tmp2
    linestr=write_echo_file(fid_inp);
    tmp=str2num(linestr);
    linestr=write_echo_file(fid_inp);
    tmp2=str2num(linestr); 
    for i=1:nnel
        if(i<=7)
           nodes(iel,i)=tmp(i+1);
        else
           nodes(iel,i)=tmp2(i-7);
        end
    end
end
%---------------------------------------------------
% determine the boundary nodes
%---------------------------------------------------
[bound]=find_node(gcoord,B);
%-----------------------------------
% input data for boundary conditions
%-----------------------------------
[bcdof,bcval]=In_bcdof(bound,ndof);
%------------------------------
% calculate system force vector
%------------------------------
[ff]=cal_ext_load(gcoord,nodes,ndof,L,p0);
%------------------------------
% input the constitutive matrix
%------------------------------
[matmtx]=fematiso(4,emodule,poisson);
%----------------------------------
% calculate system stiffness matrix 
%----------------------------------
[kk]=cal_stiff_t10(gcoord,nodes,ndof,matmtx);
%--------------------------
% apply boundary conditions
%--------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);
%---------------------------
% solve the matrix equations
%---------------------------
[LL,UU]=lu(kk);
utemp=LL\ff;
disp=UU\utemp;

solve_time=cputime-first_time
%-------------------------------------------------------------
% data post process
%-------------------------------------------------------------
x=[];
y=[];
z=[];
for i=1:nnode
    x(i)=disp(3*i-2);
    y(i)=disp(3*i-1);
    z(i)=disp(3*i);
    M(i)=sqrt(x(i).^2+y(i).^2+z(i).^2);
end
fid = fopen('Results.plt','w');
fprintf(fid,'TITLE="erjiesimianti"\n');
fprintf(fid,'VARIABLES="X""Y""Z""U""V""W""M"\n');
maxnode=max([max(nodes(:,1)),max(nodes(:,2)),max(nodes(:,3)),max(nodes(:,4))]);
fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=TETRAHEDRON,F=FEPOINT\n',maxnode,nel);
for count = 1:maxnode
    fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',gcoord(count,1),gcoord(count,2),gcoord(count,3),x(count),y(count),z(count),M(count));
end

for i=1:nel
    fprintf(fid,'%8d %8d %8d %8d\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
end

fclose(fid);