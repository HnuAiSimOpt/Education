function [Kglobal,Fglobal]=Plain_stiff(Kglobal,Fglobal,node,element,D,h,bodyforce)
% 进行系统分析，得到整体刚度矩阵

Nelem=size(element,1);

for iel=1:Nelem
    nod=element(iel,:); xx=find(nod>0);nod=nod(xx);    % 剔除零结点 
    sctr=zeros(2*length(nod),1) ;                               
    sctr(1:2:end)=2*nod-1; sctr(2:2:end)=2*nod;              % 单元的自由度
    coords=node(nod,:);                                       % 结点坐标
    center=mean(coords,1);
    
    Ke=Element_stiff_tri2D(D(center),coords,h);
    
    % 体积力
    body=bodyforce(center);area=0.5*det([ones(3,1),coords]); 
    Fe=zeros(6,1);
    Fe(1:2:6)=h*area*body(1)*[1;1;1]/3;
    Fe(2:2:6)=h*area*body(2)*[1;1;1]/3;
    
    % 单元组装
    Kglobal(sctr,sctr)=Kglobal(sctr,sctr)+Ke;
    Fglobal(sctr)=Fglobal(sctr)+Fe;
end