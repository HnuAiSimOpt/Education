function [point2,weight2]=feglqd2(nglx,ngly)

%-------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients
%     of Gauss-Legendre quadrature for two-dimensional integration
%
%  Synopsis:
%     [point2,weight2]=feglqd2(nglx,ngly) 
%
%  Variable Description:
%     nglx - number of integration points in the x-axis
%     ngly - number of integration points in the y-axis
%     point2 - vector containing integration points   
%     weight2 - vector containing weighting coefficients 
%-------------------------------------------------------------------

%  determine the largest one between nglx and ngly

   if nglx > ngly
      ngl=nglx;
   else
      ngl=ngly;
   end

%  initialization

   point2=zeros(ngl,2);
   weight2=zeros(ngl,2);

%  find corresponding integration points and weights

 [pointx,weightx]=feglqd1(nglx);     % quadrature rule for x-axis
 [pointy,weighty]=feglqd1(ngly);     % quadrature rule for y-axis

%  quadrature for two-dimension

 for intx=1:nglx                     % quadrature in x-axis
   point2(intx,1)=pointx(intx);
   weight2(intx,1)=weightx(intx);
 end

 for inty=1:ngly                     % quadrature in y-axis
   point2(inty,2)=pointy(inty);
   weight2(inty,2)=weighty(inty);
 end
  
