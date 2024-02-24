function [keG,eDofs] = transform_stiff(direcVec,keLoc2,i,eNodei,nP)
%%%%%% 变换局部刚度矩阵：keLoc2
%%%%%% 到全局坐标系
%%%%%% 定义元素并正确分配元素刚度
%%%%%% 总刚度矩阵

H0 = zeros(3,3);
H = [direcVec(i,1:3);direcVec(i,4:6);direcVec(i,7:9)];
T = [H, H0,H0,H0;
     H0,H, H0,H0;
     H0,H0,H, H0;
     H0,H0,H0,H];

keG = T'*keLoc2*T;

eDofs = [eNodei; eNodei+nP; eNodei+2*nP;...
    eNodei+3*nP; eNodei+4*nP; eNodei+5*nP];
eDofs = [eDofs(:,1);eDofs(:,2)];
end

