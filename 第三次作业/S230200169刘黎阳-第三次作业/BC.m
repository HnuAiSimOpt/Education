<<<<<<< HEAD
%% ����߽��������޸�����նȺ��غ����� 
function [KK,RR]=BC(Node_info,K,R)
num=find(Node_info(:,2)<1e-9);
KK=K;
RR=R;
for i=1:size(num,1)
 r=num(i);
 KK(2*r-1,:)=0;
 KK(:,2*r-1)=0;
 KK(2*r-1,2*r-1)=1;
 KK(2*r,:)=0;
 KK(:,2*r)=0;
 KK(2*r,2*r)=1;
 RR(2*r-1)=0;
 RR(2*r)=0;
end
=======
%% ����߽��������޸�����նȺ��غ����� 
function [KK,RR]=BC(Node_info,K,R)
num=find(Node_info(:,2)<1e-9);
KK=K;
RR=R;
for i=1:size(num,1)
 r=num(i);
 KK(2*r-1,:)=0;
 KK(:,2*r-1)=0;
 KK(2*r-1,2*r-1)=1;
 KK(2*r,:)=0;
 KK(:,2*r)=0;
 KK(2*r,2*r)=1;
 RR(2*r-1)=0;
 RR(2*r)=0;
end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
