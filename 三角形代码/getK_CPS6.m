function [kk]=getK_CPS6(gcoords,nodes,nes,Cm,nnpe,ndpn,kk)
for count = 1:nes
    ke=sparse(nnpe*ndpn,nnpe*ndpn);
    cnode = nodes(count,2:7);% 单元的节点编号
    x=zeros(length(cnode),1);
    y=zeros(length(cnode),1);
    index = feeldof(cnode,nnpe,ndpn);   % 单元节点对应的自由度
    for count2=1:length(cnode)
        x(count2)=gcoords(cnode(count2),2);
        y(count2)=gcoords(cnode(count2),3);
    end
    for GaussN=1:3
        if GaussN == 1
            L1 = 2/3; L2 = 1/6; L3 = 1/6;
        elseif GaussN == 2
            L1 = 1/6; L2 = 2/3; L3 = 1/6;
        else
            L1 = 1/6; L2 = 1/6; L3 = 2/3;
        end
        xy = [1,x(1),y(1);
              1,x(2),y(2);
              1,x(3),y(3);];
        A = det(xy)/2;
        a1 = x(2)*y(3)-x(3)*y(2);
        b1 = -y(3)+y(2);
        c1 = x(3)-x(2);
        a2 = -x(1)*y(3)+x(3)*y(1);
        b2 = y(3)-y(1);
        c2 = -x(3)+x(1);
        a3 = x(1)*y(2)-x(2)*y(1);
        b3 = -y(2)+y(1);
        c3 = x(2)-x(1);
Bmatrix = [b1*(4*L1-1)/2/A,0,b2*(4*L2-1)/2/A,0,b3*(4*L3-1)/2/A,0,2*(L2*b1+L1*b2)/A,0,2*(L3*b2+L2*b3)/A,0,2*(L3*b1+L1*b3)/A,0;
           0,c1*(4*L1-1)/2/A,0,c2*(4*L2-1)/2/A,0,c3*(4*L3-1)/2/A,0,2*(L2*c1+L1*c2)/A,0,2*(L3*c2+L2*c3)/A,0,2*(L3*c1+L1*c3)/A;
         c1*(4*L1-1)/2/A,b1*(4*L1-1)/2/A,c2*(4*L2-1)/2/A,b2*(4*L2-1)/2/A,c3*(4*L3-1)/2/A,b3*(4*L3-1)/2/A,2*(L2*c1+L1*c2)/A,2*(L2*b1+L1*b2)/A,2*(L3*c2+L2*c3)/A,2*(L3*b2+L2*b3)/A,2*(L3*c1+L1*c3)/A,2*(L3*b1+L1*b3)/A;];
     ke = ke+Bmatrix'*Cm*Bmatrix*A/3;
    end
    kk = assemblekk(kk,ke,index);
end
return