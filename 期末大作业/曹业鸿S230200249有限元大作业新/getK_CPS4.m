function [kk]=getK_CPS4(gcoord,nodes,nes,Cm,nnpe,ndpn,kk)
nglx = 2;
ngly = 2;
[PA,WA] = feglq2D(nglx,ngly);
for count=1:nes
    cnode = nodes(count,2:5);
    x = zeros(nnpe,1);
    y = zeros(nnpe,1);
    for count2=1:length(cnode)
        x(count2)=gcoord(cnode(count2),2);
        y(count2)=gcoord(cnode(count2),3);
    end
    ke=zeros(nnpe*ndpn,nnpe*ndpn);
    for countx = 1:nglx
        xi = PA(countx,1);
        wxi = WA(countx,1);
        for county = 1:ngly
            eta = PA(county,2);
            weta = WA(county,2);
            [~,dhdxi,dhdeta] = feisoq4(xi,eta);
            Jacobi = 0.25*[ eta-1, 1-eta, 1+eta, -1-eta;
                            xi-1,  -1-xi, 1+xi,  1-xi    ]*[x,y];
            detJ = det(Jacobi);
            invJ = inv(Jacobi);
            [dhdx,dhdy] = federiv2(nnpe,dhdxi,dhdeta,invJ);
            Bm = [ dhdx(1),  0,    dhdx(2), 0,     dhdx(3),  0,     dhdx(4),  0 ;
                     0,    dhdy(1),  0,    dhdy(2),  0,    dhdy(3),   0,    dhdy(4);
                   dhdy(1),dhdx(1),dhdy(2),dhdx(2),dhdy(3),dhdx(3), dhdy(4),dhdx(4);];
            ke = ke+Bm'*Cm*Bm*wxi*weta*detJ;
        end
    end
    index=feeldof(cnode,nnpe,ndpn);
    kk=assemblekk(kk,ke,index);
end
return