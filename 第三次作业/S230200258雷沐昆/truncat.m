function index=truncat(master,slave)
% 从一个大集合中消除某些指标

index=[];
N=length(master);
for i=1:N
    if ~ismember(master(i),slave)
        index=[index;i];
    end
end
