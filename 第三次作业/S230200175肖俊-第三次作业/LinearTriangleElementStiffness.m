function y = LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p)

A=(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj))/2;%���㵥Ԫ���

betai=yj-ym; %���㵥ԪӦ������е�Ԫ��

betaj=ym-yi;%���㵥ԪӦ������е�Ԫ��

betam=yi-yj;%���㵥ԪӦ������е�Ԫ��

gammai=xm-xj;%���㵥ԪӦ������е�Ԫ��

gammaj=xi-xm;%���㵥ԪӦ������е�Ԫ��

gammam=xj-xi;%���㵥ԪӦ������е�Ԫ��

B=[betai 0 betaj 0 betam 0;%�γɵ�ԪӦ�����B

    0 gammai 0 gammaj 0 gammam;

    gammai betai gammaj betaj gammam betam]/(2*A);

if p==1%����ƽ��Ӧ�������γɵ��Ծ���D

    D=(E/(1-NU^2))*[1 NU 0;NU 1 0;0 0 (1-NU)/2];

elseif p==2%����ƽ��Ӧ�������γɵ��Ծ���D

    D=(E/(1+NU)/(1-2*NU))*[1-NU NU 0;NU 1-NU 0;0 0 (1-2*NU)/2];
end

y=t*A*B'*D*B;%�����鹦ԭ���γɵ�Ԫ�նȾ���