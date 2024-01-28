function [element,node]=mesh_grid(lengthx,lengthy,d_length)
%得到各节点的几何坐标
node_num=0;
Node=[];
for i=0:d_length:lengthx
 for j=0:d_length:lengthy
     node_num=node_num+1;
     Node=[Node;node_num,i,j];
 end
end
for i=0.9:d_length:1.1
    for j=0.4:d_length:0.6
        for n=1:231
            if ((Node(n,2)==i)&&(Node(n,3)==j))
                Node(n)=0;
            end
        end
    end
end
b=0;
for n=1:231
    if Node(n,1)~=0
        b=b+1;
        node(b,1)=b;
        node(b,2)=Node(n,2);
        node(b,3)=Node(n,3);
    end
end
%得到空心结构各节点坐标
ele_num=0;
Element=[];
I=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2];
J=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
A=[0.9,1,1.1];
B=[0.4,0.5,0.6];
for i=I
 for j=J
     if (ismember(i,A)==0)||(ismember(j,B)==0)
     ele_num=ele_num+1;
     node_local=[i-d_length,j-d_length;i,j-d_length;i-d_length,j];
     node_list=[];
%得到网格单元的三个节点坐标
 for k=1:3
     row_x=find(abs(node(:,2)-node_local(k,1))<1e-6);
     row_y=find(abs(node(:,3)-node_local(k,2))<1e-6);
     num_node=intersect(row_x,row_y);
 %找到与节点坐标对应的节点编号
     node_list=[node_list;num_node];
 end 
 [m,~]=size(node_list);
     if m<3
         node_list(3)=0;
     end
 Element(ele_num,1)=ele_num;Element(ele_num,2)=node_list(1);Element(ele_num,3)=node_list(2);Element(ele_num,4)=node_list(3);
     end
  end
end
b=1;
for i=1:191
    if Element(i,4)~=0
        element(b,1)=Element(i,1);element(b,2)=Element(i,2);element(b,3)=Element(i,3);element(b,4)=Element(i,4);
        b=b+1;
    end
end
%得到每一个网格单元包含的节点编号矩阵