%% 单元应力计算 （输入应变总矩阵、单元总信息、计算得到的位移解，输出每个单元的应力）
function [sigma_x,sigma_y,sigma_xy]=Stress(BB,Ele_info,D,u)
Num_eles=size(Ele_info,1);
sigma=[];
sigma_x=[];
sigma_y=[];
sigma_xy=[];
for i=1:Num_eles
 node_local=Ele_info(i,2:4);
 u_local=[u(2*node_local(1)-1);
 u(2*node_local(1));
 u(2*node_local(2)-1);
 u(2*node_local(2));
 u(2*node_local(3)-1);
 u(2*node_local(3))];
 sigma=D*BB(3*i-2:3*i,:)*u_local;
 sigma_x=[sigma_x;sigma(1)];
 sigma_y=[sigma_y;sigma(2)];
 sigma_xy=[sigma_xy;sigma(3)];
end

