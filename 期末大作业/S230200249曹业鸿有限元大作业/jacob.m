function [jacob2]=jacob(nnel,dhdr,dhds,xcoord,ycoord)
%�õ���˹����Ÿ��Ⱦ���
%nnel�����ֵ�����

 jacob2=zeros(2,2);

 for i=1:nnel
 jacob2(1,1)=jacob2(1,1)+dhdr(i)*xcoord(i);
 jacob2(1,2)=jacob2(1,2)+dhdr(i)*ycoord(i);
 jacob2(2,1)=jacob2(2,1)+dhds(i)*xcoord(i);
 jacob2(2,2)=jacob2(2,2)+dhds(i)*ycoord(i);
 end
