function [B]=bjuzhen_1(nelem,dd)
%º∏∫Œæÿ’ÛB
 for i=1:nelem
        B(1,2*i-1)=dd(1,i);
        B(1,2*i)=0;
        B(2,2*i-1)=0;
        B(2,2*i)=dd(2,i);
        B(3,2*i-1)=dd(2,i);
        B(3,2*i)=dd(1,i); 
 end
