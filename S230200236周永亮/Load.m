%% 节点载荷列阵
function R=Load(Node_info,F)
R=zeros(size(Node_info,1)*2,1);
row1_x=find(Node_info(:,2)==2);
row1_y=find(Node_info(:,3)==0);
num1=intersect(row1_x,row1_y);
row2_x=find(Node_info(:,2)==2);
row2_y=find(Node_info(:,3)==1);
num2=intersect(row2_x,row2_y);
R(2*num1)=-F/2;
R(2*num2)=-F/2;
end