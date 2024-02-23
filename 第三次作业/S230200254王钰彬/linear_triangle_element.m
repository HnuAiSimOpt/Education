%姓名：王钰彬 学号S230200254
%本程序主要制作了线性三角形单元
clc
clear
close all
x=[1,3,5];
y=[1,1,4];
hold all;
subplot(222)
plot(x,y,'o','color','k');
patch('faces',[1,2,3],'vertices',[x;y]','facecolor','c','edgecolor','k')
x21=x(2)-x(1);
x31=x(3)-x(1);
y21=y(2)-y(1);
y31=y(3)-y(1);
e=[0,1,0];
z=[0,0,1];
hold all;
subplot(221)
plot(e,z,'o','color','k');
patch('faces',[1,2,3],'vertices',[e;z]','facecolor','c','edgecolor','k')
xm=x(1)+x21*e+x31*z;
ym=y(1)+y21*e+y31*z;
subplot(224)
plot(xm,ym,'o','color','k');
patch('faces',[1,2,3],'vertices',[xm;ym]','facecolor','c','edgecolor','k')