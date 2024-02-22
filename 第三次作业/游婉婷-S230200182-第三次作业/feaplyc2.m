function [KK,ff]=feaplyc2(node,K,load_node_matrix)
num=find(node(:,2)<1e-9);
KK=K;
ff=load_node_matrix;
for i=1:size(num,1)
 r=num(i);
 KK(2*r-1,:)=0;
 KK(:,2*r-1)=0;
 KK(2*r-1,2*r-1)=1;
 KK(2*r,:)=0;
 KK(:,2*r)=0;
 KK(2*r,2*r)=1;
 ff(2*r-1)=0;
 ff(2*r)=0;
end