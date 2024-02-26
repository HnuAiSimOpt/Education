function [point,weight]=feglqd1(n)
% 确定用于一维积分的积分点和加权系数
 if n==1          
    point(1)=0.0;
    weight(1)=2.0;
 elseif n==2     
    point(1)=-0.577350269189626;
    point(2)=-point(1);
    weight(1)=1.0;
    weight(2)=weight(1);
 elseif n==3      
    point(1)=-0.774596669241483;
    point(2)=0.0;
    point(3)=-point(1);
    weight(1)=0.555555555555556;
    weight(2)=0.888888888888889;
    weight(3)=weight(1);
 elseif n==4       
    point(1)=-0.861136311594053;
    point(2)=-0.339981043584856;
    point(3)=-point(2);
    point(4)=-point(1);
    weight(1)=0.347854845137454;
    weight(2)=0.652145154862546;
    weight(3)=weight(2);
    weight(4)=weight(1);
else                
    point(1)=-0.906179845938664;
    point(2)=-0.538469310105683;
    point(3)=0.0;
    point(4)=-point(2);
    point(5)=-point(1);
    weight(1)=0.236926885056189;
    weight(2)=0.478628670499366;
    weight(3)=0.568888888888889;
    weight(4)=weight(2);
    weight(5)=weight(1);
end
