<<<<<<< HEAD
function [kinmtx2]=feglB(nnel,dhdx,dhdy)
%��װ���������½ڵ��B����Ӧ�����

 for i=1:nnel
 i1=(i-1)*2+1;  
 i2=i1+1;
 kinmtx2(1,i1)=dhdx(i);
 kinmtx2(2,i2)=dhdy(i);
 kinmtx2(3,i1)=dhdy(i);
 kinmtx2(3,i2)=dhdx(i);
 end
=======
function [kinmtx2]=feglB(nnel,dhdx,dhdy)
%��װ���������½ڵ��B����Ӧ�����

 for i=1:nnel
 i1=(i-1)*2+1;  
 i2=i1+1;
 kinmtx2(1,i1)=dhdx(i);
 kinmtx2(2,i2)=dhdy(i);
 kinmtx2(3,i1)=dhdy(i);
 kinmtx2(3,i2)=dhdx(i);
 end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
