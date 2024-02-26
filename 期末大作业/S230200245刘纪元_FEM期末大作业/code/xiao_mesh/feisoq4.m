function [shapeq4,dhdrq4,dhdsq4]=feisoq4(rvalue,svalue)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric four-node quadilateral shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapeq4,dhdrq4,dhdsq4]=feisoq4(rvalue,svalue)  
%
%  Variable Description:
%     shapeq4 - shape functions for four-node element
%     dhdrq4 - derivatives of the shape functions w.r.t. r
%     dhdsq4 - derivatives of the shape functions w.r.t. s
%     rvalue - r coordinate value of the selected point   
%     svalue - s coordinate value of the selected point
%
%  Notes:
%     1st node at (-1,-1), 2nd node at (1,-1)
%     3rd node at (1,1), 4th node at (-1,1)
%------------------------------------------------------------------------

% shape functions

 shapeq4(1)=0.25*(1-rvalue)*(1-svalue);
 shapeq4(2)=0.25*(1+rvalue)*(1-svalue);
 shapeq4(3)=0.25*(1+rvalue)*(1+svalue);
 shapeq4(4)=0.25*(1-rvalue)*(1+svalue);

% derivatives

 dhdrq4(1)=-0.25*(1-svalue);
 dhdrq4(2)=0.25*(1-svalue);
 dhdrq4(3)=0.25*(1+svalue);
 dhdrq4(4)=-0.25*(1+svalue);

 dhdsq4(1)=-0.25*(1-rvalue);
 dhdsq4(2)=-0.25*(1+rvalue);
 dhdsq4(3)=0.25*(1+rvalue);
 dhdsq4(4)=0.25*(1-rvalue);
