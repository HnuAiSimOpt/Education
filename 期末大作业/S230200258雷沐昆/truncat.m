function index=truncat(master,slave)
% ��һ���󼯺�������ĳЩָ��

index=[];
N=length(master);
for i=1:N
    if ~ismember(master(i),slave)
        index=[index;i];
    end
end
