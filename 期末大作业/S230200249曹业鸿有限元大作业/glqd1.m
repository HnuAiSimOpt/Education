function [point1,weight1]=glqd1(ngl)
%根据高斯积分点数确定各积分点坐标和权系数

   point1=zeros(ngl,1);
   weight1=zeros(ngl,1);

 if ngl==1          
    point1(1)=0.0;
    weight1(1)=2.0;

 elseif ngl==2       
    point1(1)=0.577350269189626;
    point1(2)=-point1(1);
    weight1(1)=1.0;
    weight1(2)=weight1(1);

 elseif ngl==3       
    point1(1)=-0.774596669241483;
    point1(2)=0.0;
    point1(3)=-point1(1);
    weight1(1)=0.555555555555556;
    weight1(2)=0.888888888888889;
    weight1(3)=weight1(1);

 elseif ngl==4      
    point1(1)=-0.861136311594053;
    point1(2)=-0.339981043584856;
    point1(3)=-point1(2);
    point1(4)=-point1(1);
    weight1(1)=0.347854845137454;
    weight1(2)=0.652145154862546;
    weight1(3)=weight1(2);
    weight1(4)=weight1(1);
 
else                 
    point1(1)=-0.906179845938664;
    point1(2)=-0.538469310105683;
    point1(3)=0.0;
    point1(4)=-point1(2);
    point1(5)=-point1(1);
    weight1(1)=0.236926885056189;
    weight1(2)=0.478628670499366;
    weight1(3)=0.568888888888889;
    weight1(4)=weight1(2);
    weight1(5)=weight1(1);

end
