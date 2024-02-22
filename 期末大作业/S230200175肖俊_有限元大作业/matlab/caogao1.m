%%
%节点坐标
x0=[];
for i=1:lx+1
    for j=1:ly+1
        x0=[x0; (i-1)*lengthx/lx  -0.5*lengthy*(1+(lx+1-i)/lx)*(1-(j+16)/ly)];
    end
end

%表示各单元节点坐标，分上三角和下三角
nodes=[];%各个单元的节点坐标
for i=1:lx
    for j=1:ly
        nodes=[nodes; (ly+1)*(i-1)+j (ly+1)*i+j (ly+1)*i+j+1 (ly+1)*(i-1)+j+1;];
    end
end

%绘制网格
figure(2)
hold on
axis off%取消对坐标轴的一切设置
axis equal%严格控制各坐标的分度使其相等
for ie=1:nel%单元数
    for j=1:nnel+1%连接成环
        j1=mod(j-1,nnel)+1;%mod取余数运算（(j-1）/nnel)
        xp(j)=x0(nodes(ie,j1),1);
        yp(j)=x0(nodes(ie,j1),2);
    end
    plot(xp,yp,'-')
end