function [kinmtx2]=fekine3d(nnel,dhdx,dhdy,dhdz)


 for i=1:nnel
 i1=(i-1)*3+1;  
 i2=i1+1;
 i3=i1+2;
 kinmtx2(1,i1)=dhdx(i);
 kinmtx2(2,i2)=dhdy(i);
 kinmtx2(3,i3)=dhdz(i);
 
 kinmtx2(4,i1)=dhdy(i);
 kinmtx2(4,i2)=dhdx(i);
 
 kinmtx2(5,i2)=dhdz(i);
 kinmtx2(5,i3)=dhdy(i);
 kinmtx2(6,i1)=dhdz(i);
 kinmtx2(6,i3)=dhdx(i);
 
 
 end