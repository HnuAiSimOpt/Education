function [stress] = solveStress(num_nodes, elemNodes,...
    nodeCoordinates, constMatrix, disp, h)

%% �˺������ڼ��㵥Ԫ�ĵ�Ӧ��
% stress ȫ��Ӧ��
% num_nodes �ڵ�����
% elemNodes ��Ԫ�ڵ���
% nodeCoordinates �ڵ�����
% thickness ��Ԫ���
% constMatrix ��������
% disp �ڵ�λ������

%% 
stress = zeros(num_nodes, 6);
stress_points = [-1 -1; 1 -1; 1 1; -1 1];
for e = 1:size(elemNodes, 1)
    nodes = elemNodes(e,:);
    index = zeros(1,length(nodes)); 
    ele_disp_m = zeros(8, 1);
    ele_disp_b = zeros(8, 1);
    ele_disp_s = zeros(12, 1);
 
    for i = 1:length(nodes)
        index(i) = find(nodeCoordinates(:,1)==nodes(i));
    end
    elemcoords = nodeCoordinates(index, 2:end);
    ele_disp_m(1:2:end) = disp(6*index-5);
    ele_disp_m(2:2:end) = disp(6*index-4);
    ele_disp_b(1:2:end) = disp(6*index-2);
    ele_disp_b(2:2:end) = disp(6*index-1);
    ele_disp_s(1:3:end) = disp(6*index-3);
    ele_disp_s(2:3:end) = disp(6*index-2);
    ele_disp_s(3:3:end) = disp(6*index-1);
    
    for g = 1:size(stress_points, 1)
        xi = stress_points(g, 1);
        eta = stress_points(g, 2);
        
        [N, dN] = shapeFun2D(xi, eta, 'Q4');
        [~, dXY] = jacobian2D(elemcoords, dN);
        [Bm, Bb] = strainMatrix2D(dXY); 
        [~, ~, Bs] = strainMatrix2D(dXY, N);
        
        stress(index(g), [1 2 6]) = constMatrix.in * Bm * ele_disp_m;
        stress(index(g), [1 2 6]) = stress(index(g), [1 2 6]) + ...
            (h * constMatrix.in * Bb * ele_disp_b)';
        stress(index(g), [4 5]) = constMatrix.out * Bs * ele_disp_s;
    end
end
end