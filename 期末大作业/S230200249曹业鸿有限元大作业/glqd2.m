function [point2,weight2]=glqd2(nglx,ngly)
%确定高斯积分点坐标和对应权系数
%nglx：x方向积分点数
%ngly：y方向积分点数
%point2：积分点坐标矩阵
%weigh2：积分点权系数矩阵

   if nglx > ngly
      ngl=nglx;
   else
      ngl=ngly;
   end

   point2=zeros(ngl,2);
   weight2=zeros(ngl,2);

 [pointx,weightx]=glqd1(nglx);    
 [pointy,weighty]=glqd1(ngly);     

 for intx=1:nglx                    
   point2(intx,1)=pointx(intx);
   weight2(intx,1)=weightx(intx);
 end

 for inty=1:ngly                   
   point2(inty,2)=pointy(inty);
   weight2(inty,2)=weighty(inty);
 end
  
