function[B]=Bc(J,nnel,dndr,dnds,dndt)
% 构建单元应变矩阵B
B=zeros(6,24);%单元应变矩阵的空矩阵
J1=inv(J);
dndx=zeros(nnel,1);
dndy=zeros(nnel,1);
dndz=zeros(nnel,1);
for i=1:nnel
    dndx(i)=J1(1,1)*dndr(i)+J1(1,2)*dnds(i)+J1(1,3)*dndt(i);%形函数N对x的偏导
    dndy(i)=J1(2,1)*dndr(i)+J1(2,2)*dnds(i)+J1(2,3)*dndt(i);%形函数N对y的偏导
    dndz(i)=J1(3,1)*dndr(i)+J1(3,2)*dnds(i)+J1(3,3)*dndt(i);%形函数N对z的偏导
    Be=[dndx(i) 0 0
        0 dndy(i) 0
        0 0 dndz(i)
        dndy(i) dndx(i) 0
        0 dndz(i) dndy(i)
        dndz(i) 0 dndx(i)];
    B(:,3*(i-1)+1:3*i)=B(:,3*(i-1)+1:3*i)+Be;
end