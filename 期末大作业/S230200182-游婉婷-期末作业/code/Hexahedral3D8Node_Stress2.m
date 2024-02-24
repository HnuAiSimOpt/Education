function stress=Hexahedral3D8Node_Stress2(E,poisson,b,d,element_solid,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8)
% 该函数计算单元中心点的应力
% 输入弹性模量E，泊松比poisson
% 输入8个节点的坐标
% 输入单元的位移列阵d
% 输出单元中心的应力stress(6X1)，应力分量为Sx,Sy,Sz,Sxy,Syz,Szx
% 定义形状函数N
node=element_solid(b,2:9);
d_local=[d(3*node(1)-2);d(3*node(1)-1);d(3*node(1));d(3*node(2)-2);d(3*node(2)-1);d(3*node(2));d(3*node(3)-2);d(3*node(3)-1);d(3*node(3));d(3*node(4)-2);d(3*node(4)-1);d(3*node(4));d(3*node(5)-2);d(3*node(5)-1);d(3*node(5));d(3*node(6)-2);d(3*node(6)-1);d(3*node(6));d(3*node(7)-2);d(3*node(7)-1);d(3*node(7));d(3*node(8)-2);d(3*node(8)-1);d(3*node(8))];
Loc=[x1 y1 z1;x2 y2 z2;x3 y3 z3;x4 y4 z4;x5 y5 z5;x6 y6 z6;x7 y7 z7;x8 y8 z8;];
nx=2;
ny=2;
nz=2;
point=Feglqd(nx,ny,nz);
for ix=1:nx
    sx=point(ix,1);                    
    for iy=1:ny
        mx=point(iy,2);    
        for iz=1:nz
            tx=point(iz,3);
            dNsnt=[-(1-mx)*(1+tx)/8,-(1-sx)*(1+tx)/8,(1-sx)*(1-mx)/8; %N1-8 对 s,n,t的导数矩阵
                  (1-mx)*(1+tx)/8,-(1+sx)*(1+tx)/8,(1+sx)*(1-mx)/8;
                  (1-mx)*(1-tx)/8,-(1+sx)*(1-tx)/8,-(1+sx)*(1-mx)/8;
                  -(1-mx)*(1-tx)/8,-(1-sx)*(1-tx)/8,-(1-sx)*(1-mx)/8;
                  -(1+mx)*(1-tx)/8,(1-sx)*(1+tx)/8,(1-sx)*(1+mx)/8;
                  (1+mx)*(1+tx)/8,(1+sx)*(1+tx)/8,(1+sx)*(1+mx)/8;
                  (1+mx)*(1-tx)/8,(1+sx)*(1-tx)/8,-(1+sx)*(1+mx)/8;
                  -(1+mx)*(1-tx)/8,(1-sx)*(1-tx)/8,-(1-sx)*(1+mx)/8;];
            dNsnt=dNsnt';
            J=dNsnt*Loc;
            dNxyz=pinv(J)*dNsnt;
            B=zeros(6,24);
            for ii=1:8 % 计算B矩阵
                Bii=[dNxyz(1,ii) 0 0;0 dNxyz(2,ii) 0;0 0 dNxyz(3,ii);dNxyz(2,ii) dNxyz(1,ii) 0;0 dNxyz(3,ii) dNxyz(2,ii);dNxyz(3,ii) 0 dNxyz(1,ii);];
                B(:,3*(ii-1)+1:3*ii)=Bii;
            end 
        end
    end
end
D=E/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];
w=D*B*d_local;
% 在单元的中心处计算应力
wcent=subs(w,{sx,tx,mx},{0,0,0});
stress=double(wcent);
