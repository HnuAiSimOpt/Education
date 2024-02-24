function [point,weight]=Feglqd(nx,ny,nz)
%  判断最后一个积分点，取最长的一个积分点
   if nx > ny
      ngl=nx;
   else
      ngl=ny;
   end
% 初始化权重、积分点
   point=zeros(2,3);
   weight=zeros(2,3);
% 找到对应的积分点和权重
 [pointx,weightx]=feglqd1(nx);     % x方向积分级数     
 [pointy,weighty]=feglqd1(ny);     % y方向积分级数
 [pointz,weightz]=feglqd1(nz);     % z方向积分级数
% 二维正交
 for intx=1:nx                     
   point(intx,1)=pointx(intx);
   weight(intx,1)=weightx(intx);
 end
 for inty=1:ny                     
   point(inty,2)=pointy(inty);
   weight(inty,2)=weighty(inty);
 end 
for intz=1:nz                    
   point(intz,3)=pointz(intz);
   weight(intz,3)=weightz(intz);
end