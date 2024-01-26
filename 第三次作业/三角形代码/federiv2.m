function [dhdx,dhdy]=federiv2(nnel,dhdxi,dhdeta,invjacob)
 for i=1:nnel
    dhdx(i)=invjacob(1,1)*dhdxi(i)+invjacob(1,2)*dhdeta(i);
    dhdy(i)=invjacob(2,1)*dhdxi(i)+invjacob(2,2)*dhdeta(i);
 end
 return
