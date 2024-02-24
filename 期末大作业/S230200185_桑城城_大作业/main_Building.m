% % clear
tol = 1e-6;

%% %%%%% 1: 材料特性
E0 = 30E9;   %%%% 弹性模量N/m^2
pois = 0.2;   %%%% 泊松比

G0 = E0/2/(1+pois);%%%% 剪切模量


%% %%%%% 2: 从ANSYS导入网格

id = ~isnan(NLIST(:,1));
NLIST = NLIST(id,1:4);

id = ~isnan(ELIST(:,1));
ELIST = ELIST(id,1:8);

%%%%% 重新命名节点坐标和节点
xyz = NLIST(:,2:4);
nP = size(xyz,1);

nE = size(ELIST,1);
eNodes = ELIST(:,7:8);

%% %%%%% 3: 几何特性
%%%%% 强梁
A1 = 0.16000;
Iy1 = 0.21333E-02;
Iz1 = 0.21333E-02;
kt1 = 0.36513E-02;
ks1 = 0.84211;

%%%%%% 法向梁
A2 = 0.40000E-01;
Iy2 = 0.13333E-03;
Iz2 = 0.13333E-03;
kt2 = 0.22821E-03;
ks2 = 0.84211;


A = A1*ones(nE,1);
Iy = Iy1*ones(nE,1);
Iz = Iz1*ones(nE,1);
kt = kt1*ones(nE,1);
ks = ks1*ones(nE,1);

z_floor = (0:4:60)';%%%% 楼层的z坐标
smallB = [];
for i = 1:nE
    eNodesi = eNodes(i,:);
    xyz_ei = xyz(eNodesi,:);
    for ii = 1:length(z_floor)
        if abs(sum(xyz_ei(:,3))/2 - z_floor(ii,1)) <= tol 
            A(i,1) = A2;
            Iy(i,1) = Iy2;
            Iz(i,1) = Iz2;
            kt(i,1) = kt2;
            ks(i,1) = ks2;
            smallB = [smallB;xyz_ei];
        end
    end
end


E = E0*ones(nE,1);
G = G0*ones(nE,1);


%% %%%%%% 4: 定义单位向量
refP = [-99999999.035,212.515,0];  %%%% 参考点

[direcVec] = directVectors(nE,eNodes,xyz,refP);

%% %%%%%% 5: 刚度矩阵计算 K
nDof = 6*nP;%%%% 每个节点有6个自由度=>总自由度数=6*nP
[K] = funct_Stiff_Beam3D(nP,nE,eNodes,xyz,A,E,ks,G,Iy,Iz,kt,direcVec);
%%%%% 稀疏矩阵K = sparse(nDof,nDof)

%% %%% 6: 边界条件
fixedPs = find(abs(xyz(:,3))<= tol);

fixedDof = [fixedPs; fixedPs+nP; fixedPs+2*nP;...
    fixedPs+3*nP; fixedPs+4*nP; fixedPs+5*nP;];


%% %%% 7: 加载条件
force = zeros(nDof,1);

% % Fx = 50;
% % Fy = 0;
% % Fz = 0;
% % Mx = 0;
% % My = 0;
% % Mz = 0;
% % loadP = find(xyz(:,3)>= 32-tol & xyz(:,1)<= tol);
% % force(loadP,1) = Fx;     
% % force(loadP+1*nP,1) = Fy;
% % force(loadP+2*nP,1) = Fz;
% % force(loadP+3*nP,1) = Mx;
% % force(loadP+4*nP,1) = My;
% % force(loadP+5*nP,1) = Mz;
sigma = 50; %%%% 50 N/m^2
area = 8;
Perimeter = 2*(2+4);

fi = sigma*area/Perimeter; %%% N/m

xyzE = 0.5*(xyz(eNodes(:,1),:) + xyz(eNodes(:,2),:));
loadE = find(xyzE(:,3)>= 32-tol & xyzE(:,1)<= tol);
loadP = [];

for i = 1:length(loadE)
    ei = loadE(i,1); 
    eNodei = eNodes(ei,:); %%%%% eNodei = [node 1st, node 2nd]
    Le = sqrt((xyz(eNodei(2),1)-xyz(eNodei(1),1))^2 + ...
        (xyz(eNodei(2),2)-xyz(eNodei(1),2))^2 + ...
        (xyz(eNodei(2),3)-xyz(eNodei(1),3))^2);
    
    F_node = fi*Le/2;
    force(eNodei,1) = force(eNodei,1) + F_node;
    
    loadP = [loadP;eNodei'];
end


%% %%%% 8: 结果
disp = solution(nDof,fixedDof,K,force);



%% %%%% 9: 后处理
u_dof = (1:nP)';
v_dof = u_dof + nP;
w_dof = u_dof + 2*nP;
rx_dof = u_dof + 3*nP;
ry_dof = u_dof + 4*nP;
rz_dof = u_dof + 5*nP;

scale = 5e3;
xyz_new = xyz + scale*[disp(u_dof,1),disp(v_dof,1),disp(w_dof,1)];


%%%%% ANSYS中的变形
id = ~isnan(dispANSYS(:,1));
dispANSYS = dispANSYS(id,2:4);

id = ~isnan(rotANSYS(:,1));
rotANSYS = rotANSYS(id,2:4);

xyz_new_ANS = xyz + scale*dispANSYS;


vari = disp(u_dof,1);
vari_ANS = rotANSYS(:,2);

min_vari = min(min(vari),min(vari_ANS));
max_vari = max(max(vari),max(vari_ANS));

%%%%% FEM代码绘制结果
figure
scatterplot1(xyz_new,vari,min_vari,max_vari,10);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
view(-45,20)
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')
axis equal


%%%%% 绘制ANSYS的结果
figure
scatterplot1(xyz_new_ANS,vari_ANS,min_vari,max_vari,10);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
view(-45,20)
set(gca,'FontSize',16);
set(gca, 'FontName', 'Times New Roman')
axis equal
