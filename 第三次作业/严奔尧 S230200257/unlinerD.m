%% step 1. 定义整体刚度矩阵和节点力向量
[node_number,dummy] = size( gNode ) ;
gK = sparse( node_number * 2, node_number * 2 ) ;
gFE = sparse( node_number * 2, 1 ) ;               %整体内力向量
f = sparse( node_number * 2, 1 ) ;    
%% step 2. 计算单元刚度矩阵，并集成到整体刚度矩阵中
[element_number, dunmmy] = size( gElement ) ;
gElementStrain = zeros( element_number, 3) ;       %整体应变矩阵
gElementStress = zeros( element_number, 3) ;       %整体应力矩阵
% for ie=1:1:element_number
% k = StiffnessMatrix( ie ) ;
% AssembleStiffnessMatrix( ie, k ) ;
% end
%% step 3. 计算地面超载产生的等效节点力
[df_number,dummy] = size( gDF ) ;
for idf = 1:1:df_number
enf = EquivalentNodeForce( gDF(idf,1), gDF(idf,2), gDF(idf,3), gDF(idf,4) )*10;
i = gElement( gDF(idf,1), 1 ) ;
j = gElement( gDF(idf,1), 2 ) ;
m = gElement( gDF(idf,1), 3 ) ;
f( (i-1)*2+1 : (i-1)*2+2 ) = f( (i-1)*2+1 : (i-1)*2+2 ) + enf( 1:2 ) ;
f( (j-1)*2+1 : (j-1)*2+2 ) = f( (j-1)*2+1 : (j-1)*2+2 ) + enf( 3:4 ) ;
f( (m-1)*2+1 : (m-1)*2+2 ) = f( (m-1)*2+1 : (m-1)*2+2 ) + enf( 5:6 ) ;
end         
%% step 4 切线刚度法迭代
gDelta1=zeros(node_number * 2,1); %取初值delta0=0    
gDelta=zeros(node_number * 2,1); %取初值delta0=0    
js=0;
while true
gK=zeros( node_number * 2, node_number * 2 );
gFE=zeros(node_number * 2,1);
for ie=1:1:element_number        
delta = NodeDe( ie,gDelta1);
eps = MatrixB( ie ) * delta; %公式2求epsilon0                      
sigma0 = unlinerD(ie,eps)* (eps-gElementStrain(ie,:)')+gElementStress(ie,:)';%g公式3求sigmma0
kt = zeros( 6, 6 ) ;
h  = gMaterial( gElement(ie, 4), 3 ) ;
B = MatrixB( ie );
area = ElementArea( ie );
kt = transpose(B)*unlinerD(ie,eps)*B*h*abs(area) ;    
AssembleStiffnessMatrix( ie, kt ) ;
for ii = 1:1:3
gElementStrain(ie,ii) = eps(ii);
gElementStress(ie,ii) = sigma0(ii);
end
FE = elementforce( ie ,gElementStress) ;
AssembleFE( ie, FE ) ; 
end;
Phi=( f-gFE );
%% step 5. 处理约束条件，修改刚度矩阵和节点力向量。采用乘大数法
[bc_number,dummy] = size( gBC1 ) ;
for ibc=1:1:bc_number
n = gBC1(ibc, 1 ) ;
d = gBC1(ibc, 2 ) ;
m = (n-1)*2 + d ;
Phi(m) = gBC1(ibc, 3)* gK(m,m) * 1e20 ;
gK(m,m) = gK(m,m) * 1e20 ;
end
dd = gK \ Phi;  
gDelta=dd+gDelta;
conv=norm(Phi)/norm(f)
if js>100 || conv<1e-4              %判断是否达到精度要求
break
else
gDelta1=gDelta;     
end
end
%% step 6. 计算节点应力(采用绕节点加权平均)
gNodeStress = zeros( node_number, 6 ) ;       
for i=1:node_number                         
S = zeros( 1, 3 ) ;                         
A = 0 ;
for ie=1:1:element_number
for k=1:1:3
if i == gElement( ie, k ) 
area= ElementArea( ie ) ;
S = S + gElementStress(ie,1:3 ) * area ;
A = A + area ;
break ;
end
end
end
gNodeStress(i,1:3) = S / A ;
gNodeStress(i,6) = 0.5*sqrt( (gNodeStress(i,1)-gNodeStress(i,2))^2 + 4*gNodeStress(i,3)^2 ) ;
gNodeStress(i,4) = 0.5*(gNodeStress(i,1)+gNodeStress(i,2)) + gNodeStress(i,6) ;
gNodeStress(i,5) = 0.5*(gNodeStress(i,1)+gNodeStress(i,2)) - gNodeStress(i,6) ;
end