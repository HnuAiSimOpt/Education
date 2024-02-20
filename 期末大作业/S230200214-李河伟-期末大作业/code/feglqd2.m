function [point2,weight2]=feglqd2(nglx,ngly,nglz)
%  determine the largest one between nglx and ngly
%  �ж����һ�����ֵ㣬ȡ���һ�����ֵ�

   if nglx > ngly
      ngl=nglx;
   else
      ngl=ngly;
   end

%  initialization

   point2=zeros(2,3);
   weight2=zeros(2,3);

%  find corresponding integration points and weights   �ҵ���Ӧ�Ļ��ֵ��Ȩ��

 [pointx,weightx]=feglqd1(nglx);     % quadrature rule for x-axis      �ֱ�����˸��Ȳνڵ�����꼰Ȩ��
 [pointy,weighty]=feglqd1(ngly);     % quadrature rule for y-axis
 [pointz,weightz]=feglqd1(nglz);     % quadrature rule for z-axis


%  quadrature for two-dimension

 for intx=1:nglx                     % quadrature in x-axis
   point2(intx,1)=pointx(intx);
   weight2(intx,1)=weightx(intx);
 end

 for inty=1:ngly                     % quadrature in y-axis
   point2(inty,2)=pointy(inty);
   weight2(inty,2)=weighty(inty);
 end
  
for intz=1:nglz                     % quadrature in y-axis
   point2(intz,3)=pointz(intz);
   weight2(intz,3)=weightz(intz);
end
end