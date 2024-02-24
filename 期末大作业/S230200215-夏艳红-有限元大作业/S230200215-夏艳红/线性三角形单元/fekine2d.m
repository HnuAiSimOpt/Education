function [kinmtx2]=fekine2d(nnel,dhdx,dhdy)
 for i=1:nnel
       i1=(i-1)*2+1;  
       i2=i1+1;
       kinmtx2(1,i1)=dhdx(i);
       kinmtx2(2,i2)=dhdy(i);
       kinmtx2(3,i1)=dhdy(i);
       kinmtx2(3,i2)=dhdx(i);
 end
