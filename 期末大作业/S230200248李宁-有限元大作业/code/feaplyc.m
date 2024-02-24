function [kk,ff]=feaplyc(kk,ff,bcdof,bcval)

%----------------------------------------------------------
% 解方程[kk]{x}={ff}
%  kk 总刚度矩阵
%  ff 节点力
%  bcdof 约束的节点自由度
%  bcval 约束的节点自由度值
%

 
 n=length(bcdof);
 sdof=size(kk);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       kk(c,j)=0;
    end

    kk(c,c)=1;
    ff(c)=bcval(i);
 end

