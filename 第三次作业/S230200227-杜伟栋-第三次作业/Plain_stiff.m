function [Kglobal,Fglobal]=Plain_stiff(Kglobal,Fglobal,node,element,D,h,bodyforce)
% ����ϵͳ�������õ�����նȾ���

Nelem=size(element,1);

for iel=1:Nelem
    nod=element(iel,:); xx=find(nod>0);nod=nod(xx);    % �޳����� 
    sctr=zeros(2*length(nod),1) ;                               
    sctr(1:2:end)=2*nod-1; sctr(2:2:end)=2*nod;              % ��Ԫ�����ɶ�
    coords=node(nod,:);                                       % �������
    center=mean(coords,1);
    
    Ke=Element_stiff_tri2D(D(center),coords,h);
    
    % �����
    body=bodyforce(center);area=0.5*det([ones(3,1),coords]); 
    Fe=zeros(6,1);
    Fe(1:2:6)=h*area*body(1)*[1;1;1]/3;
    Fe(2:2:6)=h*area*body(2)*[1;1;1]/3;
    
    % ��Ԫ��װ
    Kglobal(sctr,sctr)=Kglobal(sctr,sctr)+Ke;
    Fglobal(sctr)=Fglobal(sctr)+Fe;
end