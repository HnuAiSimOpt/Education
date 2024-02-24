function [xyz,Vr,dr] = BeamGenPoints(startPoint,endPoint,nPointGen)
%%%% nPointGen = 需要生成的材质点的数量
relaVec = endPoint - startPoint;
length = (relaVec(1,1)^2+relaVec(1,2)^2+relaVec(1,3)^2)^(1/2);
Vr = 1/length*relaVec;
x = (linspace(startPoint(1,1),endPoint(1,1),nPointGen))';
y = (linspace(startPoint(1,2),endPoint(1,2),nPointGen))';
z = (linspace(startPoint(1,3),endPoint(1,3),nPointGen))';
xyz = [x,y,z];
dr = ((x(1,1)-x(2,1))^2 + (y(1,1)-y(2,1))^2 + (z(1,1)- z(2,1))^2)^(1/2);
end

