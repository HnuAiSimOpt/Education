<<<<<<< HEAD
function [dhdx,dhdy]=glb(nnel,dhdr,dhds,invjacob)
%�������������½ڵ��κ�������ƫ����

 for i=1:nnel
 dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i);
 dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i);
 end
=======
function [dhdx,dhdy]=glb(nnel,dhdr,dhds,invjacob)
%�������������½ڵ��κ�������ƫ����

 for i=1:nnel
 dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i);
 dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i);
 end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
