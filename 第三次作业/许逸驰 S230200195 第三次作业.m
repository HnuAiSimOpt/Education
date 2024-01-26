% ����Ԫ��������ҵ
% ���ݳ� S230200195
% ����������һ����Ϊ4m����1m�����0.1m����������������ģ��Ϊ200GPa�����ɱ�Ϊ0.3��
% ��˹̶����ܵ���СΪ100N/m^2��������ֱ���µľ����غɡ�
% ��ƽ�������ε�Ԫ��������Ԫ������

% ���ָ��
clear all;
close all;
clc;

% ���ϲ���
E = 200e9;          % ����ģ����Pa��
nu = 0.3;           % ���ɱ�

% ���β���
L = 4;              % ���ȣ�m��
H = 1;              % �߶ȣ�m��
T = 0.1;            % ��ȣ�m��

% �غɲ���
q = 100;            % �����غɣ�N/m^2��

% ��Ԫ����
n = 20;             % ��Ԫ����
numNodes = n + 1;   % �ڵ�����

% �����ڵ��������
x = linspace(0, L, numNodes)';
y = zeros(numNodes, 1);
z = zeros(numNodes, 1);

% ������Ԫ���Ӿ���
connectivity = [(1:numNodes-1)', (2:numNodes)'];

% ����ȫ�ָնȾ�����غ�����
K = zeros(numNodes);
F = zeros(numNodes, 1);

% ���㵥Ԫ�նȾ�����غ�����
for i = 1:n
    node1 = connectivity(i, 1);
    node2 = connectivity(i, 2);
    
    % ��Ԫ���ȺͽǶ�
    dx = x(node2) - x(node1);
    dy = y(node2) - y(node1);
    dz = z(node2) - z(node1);
    L_e = sqrt(dx^2 + dy^2 + dz^2);
    cos_theta = dx / L_e;
    sin_theta = dy / L_e;
    tan_theta = sin_theta / cos_theta;
    
    % ��Ԫ�նȾ���
    ke = (E * T / L_e) * [1, -1; -1, 1];
    
    % ��Ԫ�غ�����
    fe = (q * T * L_e / 2) * [1; 1];
    
    % ת��Ϊȫ������ϵ
    Te = [cos_theta, sin_theta; -sin_theta, cos_theta];
    Ke = Te' * ke * Te;
    Fe = Te' * fe;
    
    % ��װȫ�ָնȾ�����غ�����
    K(node1:node2, node1:node2) = K(node1:node2, node1:node2) + Ke;
    F(node1:node2) = F(node1:node2) + Fe;
end

% �߽���������˹̶���
fixedNode = 1;
K(fixedNode, :) = 0;
K(fixedNode, fixedNode) = 1;
F(fixedNode) = 0;

% �ⷽ�̣����λ������
u = K \ F;

% ����Ӧ����Ӧ��
strain = zeros(n, 1);
stress = zeros(n, 1);
for i = 1:n
    node1 = connectivity(i, 1);
    node2 = connectivity(i, 2);
    
    % ��Ԫ���ȺͽǶ�
    dx = x(node2) - x(node1);
    dy = y(node2) - y(node1);
    dz = z(node2) - z(node1);
    L_e = sqrt(dx^2 + dy^2 + dz^2);
    cos_theta = dx / L_e;
    sin_theta = dy / L_e;
    
    % ��Ԫλ������
    ue = [u(node1); u(node2)];
    
    % ��ԪӦ��
    epsilon_e = (1 / L_e) * [-cos_theta, cos_theta] * ue;
    
    % ��ԪӦ��
    sigma_e = E * epsilon_e;
    
    % �洢Ӧ���Ӧ��
    strain(i) = epsilon_e;
    stress(i) = sigma_e;
end

% ����λ��ͼ��
figure;
plot(x, u, 'b-o');
xlabel('x��m��');
ylabel('λ�ƣ�m��');
title('��������λ�Ʒֲ�');

% ����Ӧ��ͼ��
figure;
plot(x(1:end-1), strain, 'r-o');
xlabel('x��m��');
ylabel('Ӧ��');
title('��������Ӧ��ֲ�');

% ����Ӧ��ͼ��
figure;
plot(x(1:end-1), stress, 'g-o');
xlabel('x��m��');
ylabel('Ӧ����Pa��');
title('��������Ӧ���ֲ�');
