function JM4=JMV_9to4(JMV)  %网格细化
load msh_sl
k=0;

%% 二次九结点单元拆分为4个线性单元
for i=1:length(JMV(:,1))
    % 1256
    k=k+1;
    JM4(k,:)=[JMV(i,1),JMV(i,2),JMV(i,5),JMV(i,4)];
  
    % 2365
    k=k+1;
    JM4(k,:)=[JMV(i,2),JMV(i,3),JMV(i,6),JMV(i,5)];

    % 4587 

    k=k+1;
    JM4(k,:)=[JMV(i,4),JMV(i,5),JMV(i,8),JMV(i,7)];

    % 5698
    k=k+1;
    JM4(k,:)=[JMV(i,5),JMV(i,6),JMV(i,9),JMV(i,8)];
end

