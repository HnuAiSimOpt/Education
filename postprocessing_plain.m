function [strain,stress]=postprocessing_plain(disp,D,node,element)

Nelem=size(element,1);
strain=zeros(Nelem,3);
stress=zeros(Nelem,3);

for iel=1:Nelem
    nod=element(iel,:); xx=find(nod>0);nod=nod(xx);    % 剔除零结点 
    sctr=zeros(2*length(nod),1) ;                               
    sctr(1:2:end)=2*nod-1; sctr(2:2:end)=2*nod;              % 单元的自由度
    coords=node(nod,:);                                       % 结点坐标
    center=mean(coords,1);
    %三角形面积
    area=0.5*det([ones(3,1),coords]);
    index=[1,2,3,1,2];
    b=coords(index(2:4),2)-coords(index(3:5),2);
    c=-coords(index(2:4),1)+coords(index(3:5),1);
    B=zeros(3,6);
    B(1,1:2:6)=b';B(2,2:2:6)=c';
    B(3,1:2:6)=c';B(3,2:2:6)=b';
    B=B/2/area; 
    % 体积力
    strain(iel,:)=(B*disp(sctr))';
    stress(iel,:)=strain(iel,:)*D(center);
end