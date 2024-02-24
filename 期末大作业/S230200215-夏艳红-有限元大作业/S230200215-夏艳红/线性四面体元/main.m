clear
clc
format short
lengthx=2;              %length of x-axis side of problem
lengthy=2;              %length of y-axis side of proble
lengthz=2;
%h=0.2;                  %thickness of problem
elastic=2.1*10^11;            %elastic modulus
poisson=0.3;            %Poisson's ratio
iopt=4;
Nodes = load("node.txt"); %total number of nodes in system 
Load = load("load.txt")';
bdn = load("boundary.txt")';
elements = load("element.txt");
nnode=length(Nodes);%total number of nodes in system
nel=length(elements);% number of element
fload=[0,1250,0];             % the total load
%lx=16;                  % number of element in x-axis
%ly=16;                   % number of element in y-axis
%lz=16;
%nel=5*0.75*lx*ly*lz;            % number of element
nnel=8;                 %number of nodes per element
ndof=3;                 %number of dofs per node
%nnode=(lx/2+1)*(ly+1)*(lz+1)+(lx/2)*(ly+1)*(lz/2+1);    %total number of nodes in system    
sdof=nnode*ndof;        %total system dofs
edof=nnel*ndof;         %degrees of freedom per element
B=zeros(6,12,nel); 
K=zeros(sdof,sdof);

%Assemble tetrahedral units
for i=1:nel
    element_nodes_index = elements(i,2:5);
    four_nodes_matrix = Nodes(element_nodes_index,2:end);
    x1=four_nodes_matrix(1,1);
    y1=four_nodes_matrix(1,2);
    z1=four_nodes_matrix(1,3);
    x2=four_nodes_matrix(2,1);
    y2=four_nodes_matrix(2,2);
    z2=four_nodes_matrix(2,3);
    x3=four_nodes_matrix(3,1);
    y3=four_nodes_matrix(3,2);
    z3=four_nodes_matrix(3,3);
    x4=four_nodes_matrix(4,1);
    y4=four_nodes_matrix(4,2);
    z4=four_nodes_matrix(4,3);
    [k, b,D]=TetrahedronElementStiffness(elastic,poisson,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4);
    B(:,:,i)= b;
    for j = 1:4
        row_ = element_nodes_index(j);
        for m = 1:4 
            col_ = element_nodes_index(m);
            K((3*row_-2):row_*3,(3*col_-2):col_*3) = ...
                K((3*row_-2):row_*3,(3*col_-2):col_*3)...
                + k((3*j-2):3*j,(3*m-2):3*m);
        end
    end
end
F = zeros(sdof,1);
nLoad = length(Load);
for i = 1 :ndof
    F((Load-1)*ndof+i,1) = fload(i)/nLoad;
end
Constrain=constrain(bdn,ndof);
K_constrain = K;
F_constrain = F;
K_constrain(Constrain,:) = [];
F_constrain(Constrain,:) = [];
K_constrain(:,Constrain) = [];
U_ = K_constrain\F_constrain;
%Return node displacement
U=nodedisplacement(U_,bdn);
%Post-processing calculation of cell node stress
snodes = zeros(nnode, 10);
xd = U(1:3:end);
yd = U(2:3:end);
zd = U(3:3:end);
for i = 1:nel
    elements_nodes_index = elements(i,2:5);
    u = zeros(12,1);
    for index = 1:4
        u((index-1)*3+1:index*3,:) = [xd(elements_nodes_index(index));...
                                    yd(elements_nodes_index(index));...
                                    zd(elements_nodes_index(index));];
    end
    sigma = (D*B(:,:,i)*u)';
    for index = 1:4
        node_index = elements_nodes_index(index);
        snodes(node_index,1:6) = snodes(node_index,1:6) + sigma;
        snodes(node_index, 10) = snodes(node_index, 10) +1;
    end
end
snodes(:,1:6) = snodes(:,1:6)./snodes(:,10);
for i = 1:nnode
    sigma = snodes(i,1:6);
    snodes(i,7:9) =TetrahedronElementPStresses(sigma);
end
figure1 = figure();
hold on;
%figure2 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),20,snodes(:,1),'fill');
%figure3 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),20,snodes(:,2),'fill');
%figure4 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),20,snodes(:,3),'fill');
tetramesh(elements(:,2:5), Nodes(:,2:4));
hold on
%figure5 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),"g.");
axis equal;
%绘制边界节点
Boundary_nodes_coor = Nodes(bdn,2:4);
scatter3(Boundary_nodes_coor(:,1),Boundary_nodes_coor(:,2),Boundary_nodes_coor(:,3),"black*");
%绘制加载节点
Load_nodes_coor = Nodes(Load,2:4);
scatter3(Load_nodes_coor(:,1),Load_nodes_coor(:,2),Load_nodes_coor(:,3),"yo");
view([45,45])
