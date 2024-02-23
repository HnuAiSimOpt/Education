function [DeleMatrix1,DeleMatrix2] = Delete(Matrix1,Matrix2,JDZYD,ZJD,JDYS)
DeleMatrix1 = Matrix1;
DeleMatrix2 = Matrix2;
count = [];%存储约束点编号，以便同时删除多行多列
for j = 1:ZJD %遍历约束点
for i = 1:JDZYD 
if JDYS(i,j) == 1
count = [count,(j-1)*JDZYD+i];%逐步的把约束点编号存储到count_K里面
end
end
end
DeleMatrix1(count,:) = [];%划掉Matrix1的行
DeleMatrix1(:,count) = [];%划掉Matrix1的列
DeleMatrix2(count,:) = [];%划掉Matrix2的行
DeleMatrix2(:,count) = [];%划掉Matrix2的列