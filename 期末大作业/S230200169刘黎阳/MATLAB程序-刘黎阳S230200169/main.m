clc
clear 

% 
load msh_sl;
niandu=1000;  % 粘度系数1000Pa*s
u1=0.01;v1=0;
u3=0;v3=0;

% u1=0.01
JBV1=[BP1',u1*ones(size(BP1))',v1*ones(size(BP1))'];
JBV3=[BP3',u3*ones(size(BP3))',v3*ones(size(BP3))'];
JBV=[JBV1;JBV3];
% P4=1000;
P4=1000;
P2=0;

% P4=0;
JBP2=[BE2,ones(size(BE2(:,1)))*P2,ones(size(BE2(:,1)))*P2];
JBP4=[BE4,ones(size(BE4(:,1)))*P4,ones(size(BE4(:,1)))*P4];
JBP=[JBP2;JBP4];

clear JBV1 JBV3 BP1 BP3 BP4
clear JBP2 JBP4 P2 P4
clear BE1 BE3 BE4
clear u1 v1 u3 v3

B1=zeros(Nd,Nz);
B2=zeros(Nd,Nz);
D11=zeros(Nz,Nz);
D12=zeros(Nz,Nz);
D21=zeros(Nz,Nz);
D22=zeros(Nz,Nz);
C1=zeros(Nz,Nd);
C2=zeros(Nz,Nd);
F1=zeros(Nz,1);
F2=zeros(Nz,1);


for i=1:E
    for ie=1:9;
        JXYe(ie,:)=JXYV(JMV(i,ie),:);
    end 
    [Be1,Be2]=Be(JXYe);
    [De11,De12,De21,De22]=De(JXYe,niandu);
    [Ce1,Ce2]=Ce(JXYe);
    for m=1:4
        for n=1:9
            B1(JMP(i,m),JMV(i,n))=B1(JMP(i,m), JMV(i,n))+Be1(m,n);
            B2(JMP(i,m),JMV(i,n))=B2(JMP(i,m),JMV(i,n))+Be2(m,n);
        end
    end

    for m=1:9
        for n=1:9
            D11(JMV(i,m),JMV(i,n))=D11(JMV(i,m),JMV(i,n))+De11(m,n);
            D12(JMV(i,m),JMV(i,n))=D12(JMV(i,m),JMV(i,n))+De12(m,n);
            D21(JMV(i,m),JMV(i,n))=D21(JMV(i,m),JMV(i,n))+De21(m,n);
            D22(JMV(i,m),JMV(i,n))=D22(JMV(i,m),JMV(i,n))+De22(m,n);
        end
    end

    for m=1:9
        for n=1:4
            C1(JMV(i,m),JMP(i,n))=C1(JMV(i,m),    JMP(i,n))+Ce1(m,n);
            C2(JMV(i,m),JMP(i,n))=C2(JMV(i,m), JMP(i,n))+Ce2(m,n);
        end
    end
end

for i=1:length(JBP(:,1))
    for ie=1:9
        JXYe(ie,:)=JXYV(JMV(JBP(i,1),ie),:);
    end
    [Fe1,Fe2]=Fe(JXYe,JBP(i,:));
    for m=1:9
        F1(JMV(JBP(i,1),m),1)=F1(JMV(JBP(i,1),m),1)+Fe1(m,1);
        F2(JMV(JBP(i,1),m),1)=F2(JMV(JBP(i,1),m),1)+Fe2(m,1);
    end
end


K=[D11 D12 -C1
   D21 D22 -C2
    B1  B2  zeros(Nd,Nd)];
B=[-F1;-F2;zeros(Nd,1)];
N_matrix=2*Nz+Nd;

for i=1: length(JBV(:,1))
    II=JBV(i,1);
    u=JBV(i,2);
    for J=1:N_matrix
        B(J)=B(J)-K(J,II)*u;
    end

    K(II,:)=zeros(1,N_matrix);
    K(:,II)=zeros(N_matrix,1);
    K(II,II)=1;
    B(II)=u;
end

for i=1: length(JBV(:,1))
    II=Nz+JBV(i,1);
    v=JBV(i,3);
    for J=1:N_matrix
        B(J)=B(J)-K(J,II)*v;
    end

    K(II,:)=zeros(1,N_matrix);
    K(:,II)=zeros(N_matrix,1);
    K(II,II)=1;
    B(II)=v;
end

clear D11 D12 D21 D22 C1 C2 B1 B2
clear F1 F2 De11 De12 De21 De22 JXYe
clear Be1 Be2 Ce1 Ce2  JBP JMP JXYP
clear Fe1 Fe2  P_element P_side P_value
clear i i_JBP ie l_cos_theta_x m_cos_theta_y r s
x=K\B;
ux_k_1=x(1:Nz);
vy_k_1=x(1+Nz:2*Nz);
p4=x(1+2*Nz:2*Nz+Nd);
p_k_1=[Pding2Pzong(p4,JMV)]';  % 压力差值

%  输出tecplot后处理结果
E=E*4
Nz=Nz
data=[JXYV,ux_k_1,vy_k_1,sqrt(ux_k_1.^2+vy_k_1.^2),p_k_1]
JMV4=JMV_9to4(JMV)

fid_out=fopen('速度压力耦合.plt','w');
fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "pinjun_v"  "p" \n');
fprintf(fid_out,'ZONE T="flow-field", N= %5d,E=%8d,ET=QUADRILATERAL, F=FEPOINT\n',Nz,E);
format long 
data(:,4)

for i=1:Nz
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6f%16.6f%16.6e',data(i,1), data(i,2) ,data(i,3),data(i,4), data(i,5) ,data(i,6));
    fprintf(fid_out,'\n');
%    W=zeros(Nz,1)
   W(i)=data(i,3)
    X(i)=data(i,4)
    Z(i)=data(i,5)
   Y(i)=data(i,6)
end
fprintf(fid_out,'\n');
 fprintf(fid_out,'%8d%8d%8d%8d\n',JMV4);

save result_of_n1
  W=  W'
   X=  X'
     Z= Z'
        Y=  Y'
       
yuntu_u(E,W,JMV4,JXYV,1)       
yuntu_v(E,X,JMV4,JXYV,1) 
yuntu_p(E,Y,JMV4,JXYV,1)

for i=1:length(BP2)

   uu(i,1)=ux_k_1(BP2(i),1);
    yy(i)=JXYV((BP2(i)),2)*1000;
end


