function [kinmtx2]=fekine2d(nnel,dhdx,dhdy)

%------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic equation between strains and displacements
%     for two-dimensional solids
%
%  Synopsis:
%     [kinmtx2]=fekine2d(nnel,dhdx,dhdy) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*2+1;  
 i2=i1+1;
 kinmtx2(1,i1)=dhdx(i);
 kinmtx2(2,i2)=dhdy(i);
 kinmtx2(3,i1)=dhdy(i);
 kinmtx2(3,i2)=dhdx(i);
 end
