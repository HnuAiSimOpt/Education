function [sigma_x,sigma_y,tao_xy]=element_stress(B,element,matmtx,d)
elements_num=size(element,1);
sigma_x=[]; 
sigma_y=[];  
tao_xy=[];   
for i=1:elements_num
 node=element(i,2:4);
 d_local=[d(2*node(1)-1);d(2*node(1));d(2*node(2)-1);d(2*node(2));d(2*node(3)-1);d(2*node(3))];
 sigma=matmtx*B(3*i-2:3*i,:)*d_local;
 sigma_x=[sigma_x;sigma(1)];  %x方向正应力
 sigma_y=[sigma_y;sigma(2)];  %y方向正应力
 tao_xy=[tao_xy;sigma(3)];    %表面切应力
end
end