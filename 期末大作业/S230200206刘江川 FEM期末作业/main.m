%% S230200206 刘江川 有限元期末作业
%% 材料特性
E=3E9;
mu=0.2;
%% 从comsol文件中获取节点和单元
    [node_coordinate,element_node]=comsol2matlab('单元2.mphtxt',4);
    n_node=size(node_coordinate,1);
    [element_number, n_ele]=size(element_node);
%% 定义约束和载荷
    node_displacement=1-(node_coordinate(:,1)==0)*ones(1,3);
    node_load_vector=zeros(n_node,3);
    %集中载荷1，位于2500,0,400处，单位N
        node_load_vector(sum(node_coordinate==[2.5,0,0.4],2)==3,3)=node_load_vector(sum(node_coordinate==[2.5,0,0.4],2)==3,3)-10E3;
    %集中载荷2，位于2500,200,400处，单位N
        node_load_vector(sum(node_coordinate==[2.5,0.2,0.4],2)==3,3)=node_load_vector(sum(node_coordinate==[2.5,0.2,0.4],2)==3,3)-10E3;
%% 单元刚度矩阵
    [element_stiffness_matrix,B,D] = TetrahedronElementStiffiness(node_coordinate,element_node,E,mu);
%% 组装刚度矩阵
    structural_stiffness_matrix=TetrahedronAssemble(element_stiffness_matrix,element_node,node_coordinate);
%% 求解
    n_d=node_displacement';
    n_d=n_d(:);
    n_f=node_load_vector';
    n_f=n_f(:);
    K=structural_stiffness_matrix;
    for i=1:n_node
        if n_d(i)==0
            K(i,:)=0;
            K(i,i)=1;
        end
    end
    d=reshape(K\n_f,3,[])'*10;
    Strain=reshape(pagemtimes(B,reshape(d(element_node',:)',3*n_ele,1,element_number)),6,[]);
    Stress=D*Strain;
%% 可视化
    vonStress=sqrt(sum((Stress(1:3,:)-Stress([2 3 1],:)).^2/2,1));
    f1=tetramesh(element_node,node_coordinate+d,vonStress);
    node_stress=zeros(n_node,1);
    for i=1:element_number
        node_stress(element_node(i,:))=node_stress(element_node(i,:))+vonStress(i);
    end
    tal=tabulate(element_node(:));
    node_stress=node_stress./tal(:,2);
    u=node_coordinate+d;
%% 计算差值
    comsolData = importoutdata("位移和力.txt", [10, Inf]);
    [mat_nd,index1] = sortrows(node_coordinate);
    [com_nd,index2] = sortrows(comsolData);
    u_err=abs(com_nd(:,4:6)/1000-d(index1,:));
    u_err(u_err<1e-5)=0;
    u_e=u_err./(com_nd(:,4:6)+eps)*1000*100;
    u_e(isinf(u_e))=0;
    error_u=mean(u_e);
    stress_e=mean(abs(com_nd(:,7)-node_stress(index1))./com_nd(:,7)*100);