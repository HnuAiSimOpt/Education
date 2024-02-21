function [dhdx,dhdy]=glb(nnel,dhdr,dhds,invjacob)
%计算物理坐标下节点形函数及其偏导数

 for i=1:nnel
 dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i);
 dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i);
 end
