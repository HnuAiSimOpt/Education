% ����Ԫ����-�Ľڵ��
% 
clear; close all;
clc;

%% �������
E = 210000;     % ����ģ��
mu = 0.3;       % ���ɱ�
P = 10;        % �غ�

Lx = 100;       % ���� 
Ly = 30;        % ���
Thickness = 1;  % ���
h = -0.5;        % ��ʾӦ����λ�ã���СΪ��-0.5*Thickness�� 0.5*Thickness��

num_x = 100;
num_y = 30;

scaleFactor = 0.1;   % ������������
%% ������������
D = constitutiveMatrix(E, mu);

%% ��������
[~, elemnodes, coords] = rectangularMesh(Lx, Ly, num_x, num_y);
% drawingMesh(elemnodes, coords);

%% ���������ɶ�
GDof = 6 * size(coords, 1);

%% ����նȾ���
K = globalStiffness2D(GDof, elemnodes, coords, Thickness, D);

%% �߽�����
% ��߹̶�
fixedNodes = find(coords(:,2)==0);
bc = zeros(1, 5*length(fixedNodes));
bc(1:5:end) = 6 * fixedNodes - 5;
bc(2:5:end) = 6 * fixedNodes - 4;
bc(3:5:end) = 6 * fixedNodes - 3;
bc(4:5:end) = 6 * fixedNodes - 2;
bc(5:5:end) = 6 * fixedNodes - 1;
% ���ܽ½ӹ̶�
% Ty = find(coords(:,2)==0 | coords(:,2)==Lx);
% Tx = find(coords(:,3)==0 | coords(:,3)==Ly);
% Tz = find(coords(:,3)==0 | coords(:,3)==Ly | coords(:,2)==0 | coords(:,2)==Lx);
% bc = [6*Ty-5; 6*Tx-4; 6*Tz-3];


%% �غ�
% load = distributeLoad(GDof, coords, elemnodes, P );

load = sparse(GDof, 1);
loadNodes = find(coords(:,2)==Lx);
% % ����
load(6*loadNodes-4) = P*Ly / num_y;
load(6*loadNodes(1)-4) = P*Ly / num_y / 2;
load(6*loadNodes(end)-4) = P*Ly / num_y / 2;
% ����
load(6*loadNodes(1)-3) = 50;
load(6*loadNodes(end)-3) = -20;
% load(6*loadNodes-3) = P*Ly / num_y / 1e5;
% load(6*loadNodes(1)-3) = P*Ly / num_y / 2e5;
% load(6*loadNodes(end)-3) = P*Ly / num_y / 2e5;

%% ���
% λ��
disp = solveDisp(K, bc, load);
deltaDisp = zeros(size(coords,1), 6);
for i = 1:6
    deltaDisp(:,i) = disp(i:6:end);
end

% Ӧ��
stress = solveStress(size(coords, 1), elemnodes, coords, D, disp, h);

%% ����
coords(:,4) = 0;
scaleFactor = max([Lx, Ly]) * 0.2 / max(max(abs(deltaDisp(:,1:3)))); %�Զ�����
% λ��
subplot(2,1,1);
postDisplacement(elemnodes, coords, deltaDisp(:,1:3), scaleFactor, 'M', 0);
% Ӧ��
subplot(2,1,2);
postStress(elemnodes, coords, scaleFactor, deltaDisp(:,1:3), stress, 'S11', h);