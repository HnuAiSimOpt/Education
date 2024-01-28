function [NeuDOF,NeuF]=Traction_processing(node,element,TraBound,Trafun,h)

Tranode=find(TraBound(node));  % 位移边界结点编号
numT=length(Tranode);
NeuDOF=[2*Tranode-1;2*Tranode];
NeuF=zeros(2*numT,1);
% 找到面力边界上所有的单元边界
numelem=size(element,1);
for iel=1:numelem
    nod=element(iel,:); xx=find(nod>0);nod=nod(xx);    % 剔除零结点 
    nod=[nod,nod(1)];
    side=[nod(1:end-1)',nod(2:end)'];
    for si=1:size(side,1)
        sid=side(si,:);
        if ismember(sid(1),Tranode) && ismember(sid(2),Tranode)
            center=mean(node(sid,:),1);
            len=norm(node(sid(2),:)-node(sid(1),:));
            tra=Trafun(center);
            x1=find(Tranode==sid(1));
            x2=find(Tranode==sid(2));
            NeuF(x1)=NeuF(x1)+tra(1)*len*h/2;
            NeuF(x2)=NeuF(x2)+tra(1)*len*h/2;
            NeuF(numT+x1)=NeuF(numT+x1)+tra(2)*len*h/2;
            NeuF(numT+x2)=NeuF(numT+x2)+tra(2)*len*h/2;
        end
    end
end
            
        
