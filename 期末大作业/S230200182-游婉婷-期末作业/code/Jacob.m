%求形函数偏导-求解B矩阵部分 
function [J,B]=Jacob(sx,mx,tx,nel,x,y,z)
 dNdr(1)=-1/8*(1-mx)*(1-tx);  % 形函数偏导
 dNdr(2)=1/8*(1-mx)*(1-tx);
 dNdr(3)=1/8*(1+mx)*(1-tx);
 dNdr(4)=-1/8*(1+mx)*(1-tx);
 dNdr(5)=-1/8*(1-mx)*(1+tx);
 dNdr(6)=1/8*(1-mx)*(1+tx);
 dNdr(7)=1/8*(1+mx)*(1+tx);
 dNdr(8)=-1/8*(1+mx)*(1+tx);
 
 dNds(1)=-1/8*(1-sx)*(1-tx);
 dNds(2)=-1/8*(1+sx)*(1-tx);
 dNds(3)=1/8*(1+sx)*(1-tx);
 dNds(4)=1/8*(1-sx)*(1-tx);
 dNds(5)=-1/8*(1-sx)*(1+tx);
 dNds(6)=-1/8*(1+sx)*(1+tx);
 dNds(7)=1/8*(1+sx)*(1+tx);
 dNds(8)=1/8*(1-sx)*(1+tx);
 
 dNdt(1)=-1/8*(1-sx)*(1-mx);
 dNdt(2)=-1/8*(1+sx)*(1-mx);
 dNdt(3)=-1/8*(1+sx)*(1+mx);
 dNdt(4)=-1/8*(1-sx)*(1+mx);
 dNdt(5)=1/8*(1-sx)*(1-mx);
 dNdt(6)=1/8*(1+sx)*(1-mx);
 dNdt(7)=1/8*(1+sx)*(1+mx);
 dNdt(8)=1/8*(1-sx)*(1+mx);
 
 J=zeros(3,3);
 for i=1:nel   % 求解雅可比矩阵
    J(1,1)=J(1,1)+dNdr(i)*x(i);
    J(1,2)=J(1,2)+dNdr(i)*y(i);
    J(1,3)=J(1,3)+dNdr(i)*z(i);
    J(2,1)=J(2,1)+dNds(i)*x(i);
    J(2,2)=J(2,2)+dNds(i)*y(i);
    J(2,3)=J(2,3)+dNds(i)*z(i);
    J(3,1)=J(3,1)+dNdt(i)*x(i);
    J(3,2)=J(3,2)+dNdt(i)*y(i);
    J(3,3)=J(3,3)+dNdt(i)*z(i);
 end  
    inv_J=inv(J);   % 求逆
       for j=1:nel
          dNdx(j)=inv_J(1,1)*dNdr(j)+inv_J(1,2)*dNds(j)+inv_J(1,3)*dNdt(j);
          dNdy(j)=inv_J(2,1)*dNdr(j)+inv_J(2,2)*dNds(j)+inv_J(2,3)*dNdt(j);
          dNdz(j)=inv_J(3,1)*dNdr(j)+inv_J(3,2)*dNds(j)+inv_J(3,3)*dNdt(j);
       end
       for j=1:nel
          i1=(j-1)*3+1;
          i2=i1+1;
          i3=i1+2;
          B(1,i1)=dNdx(j);
          B(2,i2)=dNdy(j);
          B(3,i3)=dNdz(j);
          B(4,i1)=dNdy(j);
          B(4,i2)=dNdx(j);
          B(5,i2)=dNdz(j);
          B(5,i3)=dNdy(j); 
          B(6,i1)=dNdz(j);
          B(6,i3)=dNdx(j);
        end