<<<<<<< HEAD
function [jacob2]=jacob(nnel,dhdr,dhds,xcoord,ycoord)
%得到高斯点的雅各比矩阵
%nnel：积分点总数

 jacob2=zeros(2,2);

 for i=1:nnel
 jacob2(1,1)=jacob2(1,1)+dhdr(i)*xcoord(i);
 jacob2(1,2)=jacob2(1,2)+dhdr(i)*ycoord(i);
 jacob2(2,1)=jacob2(2,1)+dhds(i)*xcoord(i);
 jacob2(2,2)=jacob2(2,2)+dhds(i)*ycoord(i);
 end
=======
function [jacob2]=jacob(nnel,dhdr,dhds,xcoord,ycoord)
%得到高斯点的雅各比矩阵
%nnel：积分点总数

 jacob2=zeros(2,2);

 for i=1:nnel
 jacob2(1,1)=jacob2(1,1)+dhdr(i)*xcoord(i);
 jacob2(1,2)=jacob2(1,2)+dhdr(i)*ycoord(i);
 jacob2(2,1)=jacob2(2,1)+dhds(i)*xcoord(i);
 jacob2(2,2)=jacob2(2,2)+dhds(i)*ycoord(i);
 end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
