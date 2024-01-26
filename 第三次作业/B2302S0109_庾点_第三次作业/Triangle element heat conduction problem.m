% ����Ԫ���������γ���
% ѧ�ţ�B2302S0109
% �������׵�
%
% �ó�����������ȴ������⣬��������������Ԫ�ؽ�����ɢ��
% �߽���������߽��¶�Ϊ0���ұ߽��¶�Ϊ100���±߽��¶�Ϊ0��
% ʩ�ӵı߽������ڽڵ�5ʩ�Ӵ�СΪ10������
%
% ʹ�õĵ�Ԫ���ͣ�����������Ԫ��

% ��������
x = [0 1 1 0.5 0.5]; % x����
y = [0 0 1 0.5 0]; % y����

% ����ڵ�͵�Ԫ
nodes = [x' y']; % �ڵ�����
elements = [1 2 4; 2 3 4; 3 5 4]; % ��Ԫ�ڵ�����

% ����߽�����
boundaryNodes = [1 2 3]; % �߽�ڵ�����
boundaryValues = [0 100 0]; % �߽�ڵ��Ӧ���¶�ֵ

% �����������
k = 0.5; % �ȵ���

% ����ʩ����
forceNode = 5; % ʩ�����Ľڵ�����
forceMagnitude = 10; % ʩ�����Ĵ�С

% ����ϡ�����
numNodes = size(nodes, 1);
A = sparse(numNodes, numNodes);
b = zeros(numNodes, 1);

% ѭ������ÿ����Ԫ
numElements = size(elements, 1);
for i = 1:numElements
    % ��ȡ��ǰ��Ԫ�Ľڵ�����
    nodeIndices = elements(i, :);
    
    % ���㵱ǰ��Ԫ�ĸնȾ���
    [kMatrix, fVector] = computeElementStiffness(nodes(nodeIndices, :), k);
    
    % ����Ԫ�նȾ�����غ�������ӵ��ܸնȾ�����غ�������
    A(nodeIndices, nodeIndices) = A(nodeIndices, nodeIndices) + kMatrix;
    b(nodeIndices) = b(nodeIndices) + fVector;
end

% ����߽�����
A(boundaryNodes, :) = 0;
A(sub2ind(size(A), boundaryNodes, boundaryNodes)) = 1;
b(boundaryNodes) = boundaryValues;

% ʩ����
b(forceNode) = b(forceNode) + forceMagnitude;

% �����Է�����
temperatures = A\b;

% ������
disp('�ڵ��¶ȣ�');
disp(temperatures);

% ���㵥Ԫ�նȾ���ĺ���
function [kMatrix, fVector] = computeElementStiffness(nodes, k)
    % ���������
    % nodes: ��Ԫ�ڵ��������ÿ�б�ʾһ���ڵ������ [x1, y1; x2, y2; x3, y3]
    % k: �����ȵ���
    
    % ���㵥Ԫ���
    A = abs((nodes(2,1)-nodes(1,1))*(nodes(3,2)-nodes(1,2)) - (nodes(3,1)-nodes(1,1))*(nodes(2,2)-nodes(1,2))) / 2;
    
    % ���㵥Ԫ�նȾ���
    kMatrix = (k / (4 * A)) * [2, -1, -1; -1, 2, -1; -1, -1, 2];
    
    % ��ʼ����Ԫ�غ�����Ϊ������
    fVector = zeros(3, 1);
end