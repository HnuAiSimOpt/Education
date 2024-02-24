function node=Solve_node(overview)
for i=7:2982
Node(i-6,:)=strsplit(overview{i,1});  % 将元胞中每个元素的内容分来
end
Node(:,1)=[];
for i=7:2982
    node(i-6,1)=str2num(Node{i-6,1});
    node(i-6,2)=str2num(Node{i-6,2});
    node(i-6,3)=str2num(Node{i-6,3});
    node(i-6,4)=str2num(Node{i-6,4});
end
for i=1:2976  % 将坐标小于0.1的值置0
    if abs(node(i,4))<0.1
        node(i,4)=0;
    end
end
for i=1:2976
    if abs(node(i,2))<0.1
        node(i,2)=0;
    end
end
for i=1:2976
    if abs(node(i,3))<0.1
        node(i,3)=0;
    end
end