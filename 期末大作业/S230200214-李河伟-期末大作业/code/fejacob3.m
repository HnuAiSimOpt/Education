function [jacob3]=fejacob3(nnel,shapeq,dhdsq,dhdnq,dhdzq,x,y,z)



 jacob3=zeros(3,3);

 for i=1:nnel
 jacob3(1,1)=jacob3(1,1)+dhdsq(i)*x(i);       
 jacob3(1,2)=jacob3(1,2)+dhdsq(i)*y(i);
 jacob3(1,3)=jacob3(1,3)+dhdsq(i)*z(i);
 jacob3(2,1)=jacob3(2,1)+dhdnq(i)*x(i);
 jacob3(2,2)=jacob3(2,2)+dhdnq(i)*y(i);
 jacob3(2,3)=jacob3(2,3)+dhdnq(i)*z(i);
 jacob3(3,1)=jacob3(3,1)+dhdzq(i)*x(i);
 jacob3(3,2)=jacob3(3,2)+dhdzq(i)*y(i);
 jacob3(3,3)=jacob3(3,3)+dhdzq(i)*z(i);
 
 end
 
 
 
 end
