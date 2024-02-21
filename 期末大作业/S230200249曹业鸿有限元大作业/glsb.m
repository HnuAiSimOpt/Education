 function [shapeq4,dhdrq4,dhdsq4]=glsb(rvalue,svalue)
%����õ���˹���ֵ���κ�����ƫ����
%shapeq4����˹���κ���
%dhdrq4����˹���κ�������x�����ƫ����
%dhdsq4����˹���κ�������y�����ƫ����
%rvalue����˹��x����
%svalue����˹��y����

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
