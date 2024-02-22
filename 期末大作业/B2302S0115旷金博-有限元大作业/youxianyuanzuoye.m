syms Xi Eta

% coordinate
XX = [0 0 1 1]; YY = [1 0 0 1];
for i = 1:4
    XX(1,i) = XX(i); XY(2,i)=YY(i);
end

% Material property; Force magnitude; geometry
v=0.25;E=210;F=10;alength=1;

%boundary conditions
%the DOF where is not restrained but unknow is set as 1.0
ux=[0 0 1 0]; uy=[0 0 1 1];
for i=1:4
    uxy(2*i-1)=ux(i);uxy(2*i)=uy(i);
end

% external force 
% the force where is unknow is set as 0.0
force=[0 0 0 0 F/sqrt(2) -F/sqrt(2) 0 0];

% shape function
N(1)=(1/4)*(1-Xi)*(1-Eta);
N(2)=(1/4)*(1+Xi)*(1-Eta);
N(3)=(1/4)*(1+Xi)*(1+Eta);
N(4)=(1/4)*(1-Xi)*(1+Eta);

% Jacobain
% compute the integration point in node 2 where Xi=0 and Eta=0;
for i =1:4
    dNdXi=diff(N,Xi); dNdEta=diff(N,Eta);
    dNdXiEta(i,1)=dNdXi(i); dNdXiEta(i,2)=dNdEta(i);
end

dNdXiEta(Xi,Eta)=dNdXiEta;

Jaco = XY*dNdXiEta(0,0);
detJ = det(Jaco)
Jinv = inv(Jaco)

dNdXY=dNdXiEta(0,0)*Jinv;

for i = 1:4
    dNdX(i) = dNdXY(i,1);dNdY(i) = dNdXY(i,2)
end

% [B] matrix (Strain relation matrix)
BMat=zeros(3,8);
for i= 1:4
    BMat(1,2*i-1) = dNdX(i);
    BMat(2,2*i) = dNdY(i);
    BMat(3,2*i-1) = dNdY(i); BMat(3,2*i) = dNdX(i);
end

% Local elastic stiffness
lam = v*E/(1+v)/(1-2*v); mu = E/2/(1+v);
CMat = [lam+2*mu lam 0;
        lam lam+2*mu 0;
        0 0 mu];

% stiffness matrix
kMat = BMat.' * CMat * BMat * detJ *2 *2

% elimation of stiffness matrix
j = 0;
for i =1:8
    if (uxy(i) == 0)
        j = j+1;
        igood(j) = i;
    end
end

for j1 = 1:j
    for j2 = 1:j
        kNew(j1,j2) = kMat(igood(j1), igood(j2));
    end
end

% elimation of force 
for j1 =1:j
    forceNew(j1) = force(igood(j1));
end

% elimation of force 
for j1=1:j
    forceNew(j1) = force(igood(j1));
end

% calculate displacement
uNew = vpa(inv(kNew)*forceNew.',4)

for j1= 1:j
    uxy(igood(j1)) = uNew(j1);
end

% calculate force 
force = vpa(kMat*uxy.',4)


% calculate strain & stress
strain = BMat*uxy.' ; stress = CMat * strain



