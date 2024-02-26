global GC T ke a F f q L Kn;
E=7e6;
A=30;
t=0.1;
G=E/2.6;
p=fopen('inputdata.txt','r');
JDS=fscanf(p,'%d',1);%%节点数
JDF=fscanf(p,'%d',1);%%每个节点自由度
JDZB=zeros(JDS,2);%%生成节点坐标矩阵
for i=1:JDS
    JDZB(i,:)=fscanf(p,'%d',2);
end 
GS=fscanf(p,'%d',1);%%杆单元总数
GBH=zeros(GS,2);%%杆单元编号矩阵
for i=1:GS
    GBH(i,:)=fscanf(p,'%d',2);
end
BS=fscanf(p,'%d',1);%%板数
BBH=zeros(BS,4);%%板单元编号数组
for i=1:BS
    BBH(i,:)=fscanf(p,'%d',4);
end 
YS=fscanf(p,'%d',1);%%约束数
YSS=zeros(YS,2);%%约束数值
for i=1:YS
    YSS(i,:)=fscanf(p,'%d',2);
end 
ZH=fscanf(p,'%d',1);%%载荷数
ZHS=zeros(ZH,2);%%载荷数值
for i=1:ZH
    ZHS(i,:)=fscanf(p,'%d',2);
end 
N=JDS*JDF;
K=zeros(N,N);
Ke2=zeros(8,8);
for i=1:GS
   GC(i)=sqrt((JDZB(GBH(i,2),1)-JDZB(GBH(i,1),1))^2+(JDZB(GBH(i,2),2)-JDZB(GBH(i,1),2))^2);
   %%杆长
   V=[JDZB(GBH(i,2),1)-JDZB(GBH(i,1),1),JDZB(GBH(i,2),2)-JDZB(GBH(i,1),2)]/GC(i);
   %%杆单位向量
   ax=V*[1;0];ay=V*[0;1];
   T(:,:,i)=[ax,ay,0,0;0,0,ax,ay];
   ke(:,:,i)=E*A/GC(i)*[1,-1;-1,1];%%杆局部单元刚度矩阵
   Ke=T(:,:,i)'*ke(:,:,i)*T(:,:,i);%%整体坐标下单元刚度矩阵
   K(2*GBH(i,1)-1:2*GBH(i,1),2*GBH(i,1)-1:2*GBH(i,1))=K(2*GBH(i,1)-1:2*GBH(i,1),2*GBH(i,1)-1:2*GBH(i,1))+Ke(1:2,1:2);
   K(2*GBH(i,1)-1:2*GBH(i,1),2*GBH(i,2)-1:2*GBH(i,2))=K(2*GBH(i,1)-1:2*GBH(i,1),2*GBH(i,2)-1:2*GBH(i,2))+Ke(1:2,3:4);
   K(2*GBH(i,2)-1:2*GBH(i,2),2*GBH(i,1)-1:2*GBH(i,1))=K(2*GBH(i,2)-1:2*GBH(i,2),2*GBH(i,1)-1:2*GBH(i,1))+Ke(3:4,1:2);
   K(2*GBH(i,2)-1:2*GBH(i,2),2*GBH(i,2)-1:2*GBH(i,2))=K(2*GBH(i,2)-1:2*GBH(i,2),2*GBH(i,2)-1:2*GBH(i,2))+Ke(3:4,3:4);
end
 for i=1:BS
    a1=JDZB(BBH(i,1),1)-JDZB(BBH(i,3),1);
    a2=JDZB(BBH(i,1),2)-JDZB(BBH(i,3),2);
    a3=JDZB(BBH(i,4),1)-JDZB(BBH(i,2),1);
    a4=JDZB(BBH(i,4),2)-JDZB(BBH(i,2),2);
    a5=-a1;
    a6=-a2;
    a7=-a3;
    a8=-a4;
    a(:,:,i)=[a1,a2,a3,a4,a5,a6,a7,a8];%%板单元a矩阵
    F(i)=1/2*sqrt((a5*a4-a3*a6)^2);%%板面积
    Ke2=a(:,:,i)'*G*t/4/F(i)*a(:,:,i);
   
    K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,1)-1:2*BBH(i,1))=K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,1)-1:2*BBH(i,1))+Ke2(1:2,1:2);
    K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,2)-1:2*BBH(i,2))=K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,2)-1:2*BBH(i,2))+Ke2(1:2,3:4);
    K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,3)-1:2*BBH(i,3))=K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,3)-1:2*BBH(i,3))+Ke2(1:2,5:6);
    K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,4)-1:2*BBH(i,4))=K(2*BBH(i,1)-1:2*BBH(i,1),2*BBH(i,4)-1:2*BBH(i,4))+Ke2(1:2,7:8);
    
    K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,1)-1:2*BBH(i,1))=K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,1)-1:2*BBH(i,1))+Ke2(3:4,1:2);
    K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,2)-1:2*BBH(i,2))=K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,2)-1:2*BBH(i,2))+Ke2(3:4,3:4);
    K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,3)-1:2*BBH(i,3))=K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,3)-1:2*BBH(i,3))+Ke2(3:4,5:6);
    K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,4)-1:2*BBH(i,4))=K(2*BBH(i,2)-1:2*BBH(i,2),2*BBH(i,4)-1:2*BBH(i,4))+Ke2(3:4,7:8);
    
    K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,1)-1:2*BBH(i,1))=K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,1)-1:2*BBH(i,1))+Ke2(5:6,1:2);
    K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,2)-1:2*BBH(i,2))=K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,2)-1:2*BBH(i,2))+Ke2(5:6,3:4);
    K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,3)-1:2*BBH(i,3))=K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,3)-1:2*BBH(i,3))+Ke2(5:6,5:6);
    K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,4)-1:2*BBH(i,4))=K(2*BBH(i,3)-1:2*BBH(i,3),2*BBH(i,4)-1:2*BBH(i,4))+Ke2(5:6,7:8);
    
    K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,1)-1:2*BBH(i,1))=K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,1)-1:2*BBH(i,1))+Ke2(7:8,1:2);
    K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,2)-1:2*BBH(i,2))=K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,2)-1:2*BBH(i,2))+Ke2(7:8,3:4);
    K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,3)-1:2*BBH(i,3))=K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,3)-1:2*BBH(i,3))+Ke2(7:8,5:6);
    K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,4)-1:2*BBH(i,4))=K(2*BBH(i,4)-1:2*BBH(i,4),2*BBH(i,4)-1:2*BBH(i,4))+Ke2(7:8,7:8);
 end 
 P=zeros(2*JDS,1);
 for i=1:ZH
     P(ZHS(i,1))=ZHS(i,2);
 end%%形成节点力矩阵
 for i=1:YS
     j=YSS(i,1);
     K(j,j)=10^40;
 end%%置大数法处理总刚度矩阵
 for i=1:YS
     j=YSS(i,1);
     P(j)=10^40*YSS(i,2);
 end%%相应处理力矩阵
 Kn=K^(-1);
 u=Kn*P;
 fprintf('各自由度位移矩阵为：\n');
 fprintf('%1.4f\n',u);%%输出位移矩阵
 fprintf('各板剪流：\n');
 ui1=zeros(8,1);
 
 for i=1:BS
     ui1=[u(2*BBH(i,1)-1);u(2*BBH(i,1));u(2*BBH(i,2)-1);u(2*BBH(i,2));u(2*BBH(i,3)-1);u(2*BBH(i,3));u(2*BBH(i,4)-1);u(2*BBH(i,4))];
    q(i)=-G*t/2/F(i)*a(:,:,i)*ui1;%%板受到的剪流
     fprintf('%3.4f\n',q(i));
 end
 f=zeros(2,1);
 L=zeros(2,1);
 for i=1:GS%%杆数
     ui2=[u(2*GBH(i,1)-1);u(2*GBH(i,1));u(2*GBH(i,2)-1);u(2*GBH(i,2))];%%杆四个节点对应的位移
   
     ui=T(:,:,i)*ui2;%%转换成局部坐标下位移
     fprintf('第“%d”杆局部坐标下位移\n',i);
     fprintf('%3.4f\n',ui);
     f=ke(:,:,i)*ui;%%杆端力
    fprintf('第“%d”杆等轴力杆力\n',i);
     fprintf('%3.4f\n',f);
     L=zeros(2,1);
             for m=1:BS
                 if((GBH(i,1)==BBH(m,1)||GBH(i,1)==BBH(m,2)||GBH(i,1)==BBH(m,3)||GBH(i,1)==BBH(m,4))&&(GBH(i,2)==BBH(m,1)||GBH(i,2)==BBH(m,2)||GBH(i,2)==BBH(m,3)||GBH(i,2)==BBH(m,4)))%%找杆对应板
                     
                     
                     if((GBH(i,1)==BBH(m,1))&&(GBH(i,2)==BBH(m,4)))%%杆对应板的腰，第一个点
                        
                            L(1)=L(1)-1/2*q(m)*GC(i);%%节点力
                            L(2)=L(2)+1/2*q(m)*GC(i);%%节点力
                    
                            
                      elseif(GBH(i,1)==BBH(m,2)&&GBH(i,2)==BBH(m,3))%%杆对应板的腰，第二个点
                            L(1)=L(1)+1/2*q(m)*GC(i);%%节点力
                            L(2)=L(2)-1/2*q(m)*GC(i);%%节点力
                      elseif(GBH(i,1)==BBH(m,3)&&GBH(i,2)==BBH(m,2))%%杆对应板的腰，第三个点
                            L(1)=L(1)-1/2*q(m)*GC(i);%%节点力
                            L(2)=L(2)+1/2*q(m)*GC(i);%%节点力
                      elseif(GBH(i,1)==BBH(m,4)&&GBH(i,2)==BBH(m,1))%%杆对应板的腰，第四个点
                            L(1)=L(1)+1/2*q(m)*GC(i);%%节点力
                            L(2)=L(2)-1/2*q(m)*GC(i);%%节点力      
                      elseif(GBH(i,1)==BBH(m,1)&&GBH(i,2)==BBH(m,2))%%杆对应板的平行边第一个点
                               
                          for j=1:GS 
                                       if ((BBH(m,3)==GBH(j,1)&&BBH(m,4)==GBH(j,2))||(BBH(m,4)==GBH(j,1)&&BBH(m,3)==GBH(j,2)))
                                         L(1)=L(1)-1/2*q(m)*GC(j);   
                                         L(2)=L(2)+1/2*q(m)*GC(j);
                                       end
                          end
                          
                        
                      elseif(GBH(i,1)==BBH(m,2)&&GBH(i,2)==BBH(m,1))%%杆对应板的平行边第二个点
                          
                          for j=1:GS 
                                       if ((BBH(m,3)==GBH(j,1)&&BBH(m,4)==GBH(j,2))||(BBH(m,4)==GBH(j,1)&&BBH(m,3)==GBH(j,2)))
                                         L(1)=L(1)+1/2*q(m)*GC(j);   
                                         L(2)=L(2)-1/2*q(m)*GC(j);
                                       end
                          end
                          
                          
                      
                      elseif(GBH(i,1)==BBH(m,3)&&GBH(i,2)==BBH(m,4))%%杆对应板的平行边第三点
                                              
                          for j=1:GS 
                                       if ((BBH(m,1)==GBH(j,1)&&BBH(m,2)==GBH(j,2))||(BBH(m,2)==GBH(j,1)&&BBH(m,1)==GBH(j,2)))
                                         L(1)=L(1)-1/2*q(m)*GC(j);   
                                         L(2)=L(2)+1/2*q(m)*GC(j);
                                       end
                          end
                          
                          
                      elseif(GBH(i,1)==BBH(m,4)&&GBH(i,2)==BBH(m,3))%%杆对应板的平行边第四个点
                                         
                          for j=1:GS 
                                       if ((BBH(m,1)==GBH(j,1)&&BBH(m,2)==GBH(j,2))||(BBH(m,2)==GBH(j,1)&&BBH(m,1)==GBH(j,2)))
                                         L(1)=L(1)+1/2*q(m)*GC(j);   
                                         L(2)=L(2)-1/2*q(m)*GC(j);
                                       end
                          end
                     end
                 end
               
             end
             L(1)=L(1)-f(1);
             L(2)=L(2)+f(2);
             
        fprintf('第“%d”杆端力\n',i);  
     fprintf('%3.4f\n',L);
 end
 
 
 
 

