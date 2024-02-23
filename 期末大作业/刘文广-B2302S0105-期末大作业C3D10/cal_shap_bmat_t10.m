function [shap,bmat]=cal_shap_bmat_t10(nnel,ndof,L,a,b,c,d,vol)
%--------------------------------------------------------------------------
% Purpose:
%     Calculate element shape function and kinematic matrix 
%
% Synopsis:
%     [shap,bmat]=cal_shap_bmat_t10(nnel,ndof,L,a,b,c,d,vol)
%
% Variable descriptions:
%     shap - element shape function 
%     bmat - element kinematic matrix 
%     nnel - number of node per element
%     ndof - number of dofs per node 
%     L - volume coordinates
%     vol - element volume
%--------------------------------------------------------------------------
edof=nnel*ndof;
shap=zeros(3,edof);
bmat=zeros(6,edof);

shap(1,1)=(2*L(1)-1)*L(1);
shap(1,4)=(2*L(2)-1)*L(2);
shap(1,7)=(2*L(3)-1)*L(3);
shap(1,10)=(2*L(4)-1)*L(4);
shap(1,13)=4*L(1)*L(2);
shap(1,16)=4*L(2)*L(3);
shap(1,19)=4*L(1)*L(3);
shap(1,22)=4*L(1)*L(4);
shap(1,25)=4*L(2)*L(4);
shap(1,28)=4*L(3)*L(4);

%---------------------------------------------

pa=1/(6*vol);
pb=2/(3*vol);

bmat(1,1)=pa*b(1)*(4*L(1)-1);
bmat(2,2)=pa*c(1)*(4*L(1)-1);
bmat(3,3)=pa*d(1)*(4*L(1)-1);

bmat(1,4)=pa*b(2)*(4*L(2)-1);
bmat(2,5)=pa*c(2)*(4*L(2)-1);
bmat(3,6)=pa*d(2)*(4*L(2)-1);

bmat(1,7)=pa*b(3)*(4*L(3)-1);
bmat(2,8)=pa*c(3)*(4*L(3)-1);
bmat(3,9)=pa*d(3)*(4*L(3)-1);

bmat(1,10)=pa*b(4)*(4*L(4)-1);
bmat(2,11)=pa*c(4)*(4*L(4)-1);
bmat(3,12)=pa*d(4)*(4*L(4)-1);

bmat(1,13)=pb*(b(2)*L(1)+b(1)*L(2));
bmat(2,14)=pb*(c(2)*L(1)+c(1)*L(2));
bmat(3,15)=pb*(d(2)*L(1)+d(1)*L(2));

bmat(1,16)=pb*(b(3)*L(2)+b(2)*L(3));
bmat(2,17)=pb*(c(3)*L(2)+c(2)*L(3));
bmat(3,18)=pb*(d(3)*L(2)+d(2)*L(3));

bmat(1,19)=pb*(b(3)*L(1)+b(1)*L(3));
bmat(2,20)=pb*(c(3)*L(1)+c(1)*L(3));
bmat(3,21)=pb*(d(3)*L(1)+d(1)*L(3));

bmat(1,22)=pb*(b(4)*L(1)+b(1)*L(4));
bmat(2,23)=pb*(c(4)*L(1)+c(1)*L(4));
bmat(3,24)=pb*(d(4)*L(1)+d(1)*L(4));

bmat(1,25)=pb*(b(4)*L(2)+b(2)*L(4));
bmat(2,26)=pb*(c(4)*L(2)+c(2)*L(4));
bmat(3,27)=pb*(d(4)*L(2)+d(2)*L(4));

bmat(1,28)=pb*(b(4)*L(3)+b(3)*L(4));
bmat(2,29)=pb*(c(4)*L(3)+c(3)*L(4));
bmat(3,30)=pb*(d(4)*L(3)+d(3)*L(4));

for i=1:nnel
    i1=(i-1)*ndof+1;
    i2=i1+1;
    i3=i2+1;
    shap(2,i2)=shap(1,i1);
    shap(3,i3)=shap(1,i1);
    
    bmat(4,i1)=bmat(2,i2);
    bmat(4,i2)=bmat(1,i1);
    bmat(5,i2)=bmat(3,i3);
    bmat(5,i3)=bmat(2,i2);
    bmat(6,i1)=bmat(3,i3);
    bmat(6,i3)=bmat(1,i1);
end

return
