function [dhdx,dhdy,dhdz]=federiv2(nnel,dhdsq,dhdnq,dhdzq,invjacob)



 for i=1:nnel
 dhdx(i)=invjacob(1,1)*dhdsq(i)+invjacob(1,2)*dhdnq(i)+invjacob(1,3)*dhdzq(i);
 dhdy(i)=invjacob(2,1)*dhdsq(i)+invjacob(2,2)*dhdnq(i)+invjacob(2,3)*dhdzq(i);
 dhdz(i)=invjacob(3,1)*dhdsq(i)+invjacob(3,2)*dhdnq(i)+invjacob(3,3)*dhdzq(i);
 end
