function[B]=Strain(J,nnel,dndr,dnds,dndt)
% 构建单元应变矩阵
B=zeros(6,24);
J1=inv(J);
dndx=zeros(nnel,1);
dndy=zeros(nnel,1);
dndz=zeros(nnel,1);
for i=1:nnel
    dndx(i)=J1(1,1)*dndr(i)+J1(1,2)*dnds(i)+J1(1,3)*dndt(i);
    dndy(i)=J1(2,1)*dndr(i)+J1(2,2)*dnds(i)+J1(2,3)*dndt(i);
    dndz(i)=J1(3,1)*dndr(i)+J1(3,2)*dnds(i)+J1(3,3)*dndt(i);
    Be=[dndx(i) 0 0
        0 dndy(i) 0
        0 0 dndz(i)
        dndy(i) dndx(i) 0
        0 dndz(i) dndy(i)
        dndz(i) 0 dndx(i)];
    B(:,3*(i-1)+1:3*i)=B(:,3*(i-1)+1:3*i)+Be;
end