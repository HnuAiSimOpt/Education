function kk = getK_S4R(gcoords,nodes,nes,Cm,Cb,Cs,t,nnpe,ndpn,kk)
for count = 1:nes
    cnode = nodes(count,2:5);
    x = zeros(nnpe,1);
    y = zeros(nnpe,1);
    z = zeros(nnpe,1);
    for count2=1:length(cnode)
        x(count2)=gcoords(cnode(count2),2);
        y(count2)=gcoords(cnode(count2),3);
        z(count2)=gcoords(cnode(count2),4);
    end
    [T3,T24] = getTmatrix(x,y,z);
    temp = T3*[x,y,z]';
    x = temp(1,:)';
    y = temp(2,:)';
    nglx = 2;
    ngly = 2;
    [PA,WA] = feglq2D(nglx,ngly);
    kem = sparse(0);
    keb = sparse(0);
    kes = sparse(0);
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
            Bm = [ dhdx(1),  0,    0,0,0,0, dhdx(2),  0,    0,0,0,0, dhdx(3),  0,    0,0,0,0, dhdx(4),  0,    0,0,0,0;
                     0,    dhdy(1),0,0,0,0,  0,     dhdy(2),0,0,0,0,    0,   dhdy(3),0,0,0,0,   0,    dhdy(4),0,0,0,0;
                   dhdy(1),dhdx(1),0,0,0,0, dhdy(2),dhdx(2),0,0,0,0, dhdy(3),dhdx(3),0,0,0,0, dhdy(4),dhdx(4),0,0,0,0;];
            kem = kem+Bm'*Cm*Bm*wxi*weta*detJ*t;
            Bb = [0,0, 0, -dhdx(1),    0    ,0,0,0, 0, -dhdx(2),    0    ,0,0,0, 0, -dhdx(3),    0    ,0,0,0, 0, -dhdx(4),    0    ,0;
                  0,0, 0,      0 ,  -dhdy(1),0,0,0, 0,      0 ,  -dhdy(2),0,0,0, 0,      0 ,  -dhdy(3),0,0,0, 0,      0 ,  -dhdy(4),0;
                  0,0, 0, -dhdy(1), -dhdx(1),0,0,0, 0, -dhdy(2), -dhdx(2),0,0,0, 0, -dhdy(3), -dhdx(3),0,0,0, 0, -dhdy(4), -dhdx(4),0,];
            keb = keb+Bb'*Cb*Bb*wxi*weta*detJ*t;
        end
    end
    nglxs = 1;
    nglys = 1;
    [PAS,WAS] = feglq2D(nglxs,nglys);
    for countx = 1:nglxs
        xi = PAS(countx,1);
        wxi = WAS(countx,1);
        for county = 1:nglys
            eta = PAS(county,2);
            weta = WAS(county,2);
            [shapeq4,dhdxi,dhdeta] = feisoq4(xi,eta);
            Jacobi = 0.25.*[ eta-1, 1-eta, 1+eta, -1-eta;
                             xi-1,  -1-xi,  1+xi,  1-xi  ]*[x,y];
            detJ = det(Jacobi);
            invJ = inv(Jacobi);
            [dhdx,dhdy] = federiv2(nnpe,dhdxi,dhdeta,invJ);
            Bs = [ 0,0,-dhdx(1),-shapeq4(1),0,0,0,0,-dhdx(2),-shapeq4(2),0,0,0,0,-dhdx(3),-shapeq4(3),0,0,0,0,-dhdx(4),-shapeq4(4),0,0;
                   0,0,-dhdy(1),0,-shapeq4(1),0,0,0,-dhdy(2),0,-shapeq4(2),0,0,0,-dhdy(3),0,-shapeq4(3),0,0,0,-dhdy(4),0,-shapeq4(4),0;];

            kes = kes+Bs'*Cs*Bs*wxi*weta*detJ*t;
        end
    end
    flag = 0;
    if length(kem(:,1))==length(kem(1,:))
        if length(keb(:,1))==length(keb(1,:))
            if length(kes(:,1))==length(kes(1,:))
                if length(keb(:,1))==length(kes(1,:))
                    flag = 1;
                end
            end
        end
    end
    if flag == 0
        fprintf('Err:The dimension of stiffness matrix is wrong!\n');
    end
    ke = kem+keb+kes;
    for c=1:4
        ke(c*6,c*6)=1;
    end
    Ke = T24'*ke*T24;
    index = feeldof(cnode,nnpe,ndpn);
    [kk] = assemblekk(kk,Ke,index);
    dispcp(count,nes);
end
return