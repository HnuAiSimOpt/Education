 function [shapeq4,dhdrq4,dhdsq4]=glsb(rvalue,svalue)
%计算得到高斯积分点的形函数及偏倒数
%shapeq4：高斯点形函数
%dhdrq4：高斯点形函数关于x方向的偏导数
%dhdsq4：高斯点形函数关于y方向的偏导数
%rvalue：高斯点x坐标
%svalue：高斯点y坐标

 shapeq4(3)=0.25*(1-rvalue)*(1-svalue);
 shapeq4(4)=0.25*(1+rvalue)*(1-svalue);
 shapeq4(1)=0.25*(1+rvalue)*(1+svalue);
 shapeq4(2)=0.25*(1-rvalue)*(1+svalue);

 dhdrq4(3)=-0.25*(1-svalue);
 dhdrq4(4)=0.25*(1-svalue);
 dhdrq4(1)=0.25*(1+svalue);
 dhdrq4(2)=-0.25*(1+svalue);

 dhdsq4(3)=-0.25*(1-rvalue);
 dhdsq4(4)=-0.25*(1+rvalue);
 dhdsq4(1)=0.25*(1+rvalue);
 dhdsq4(2)=0.25*(1-rvalue);
