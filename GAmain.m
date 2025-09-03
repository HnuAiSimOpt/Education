 clear all

data1=xlsread("D:\弟弟\只用时域时频域特征\data_x1_PCA.xlsx");
 data2=xlsread("D:\弟弟\只用时域时频域特征\data_x4_PCA.xlsx");
 data3=xlsread("D:\弟弟\只用时域时频域特征\data_x6_PCA.xlsx");
 data11=data1(2:316,2:4);
 data12=data2(2:316,2:4);
 data13=data3(2:316,2:4);

 normalized_data1 = zscore(data11);
 normalized_data2 = zscore(data12);
 normalized_data3 = zscore(data13);
 train_wine_labels=[ones(55,1); 2*ones(165,1); 3*ones(95,1);ones(55,1); 2*ones(165,1); 3*ones(95,1)];
train_wine= [normalized_data1; normalized_data2;];
ga_option.maxgen = 100;
ga_option.sizepop = 10; 
ga_option.cbound = [0,100];
ga_option.gbound = [0,100];
ga_option.v = 5;
ga_option.ggap = 0.9;
[bestacc,bestc,bestg] = gaSVMcgForClass(train_wine_labels,train_wine,ga_option);
str = sprintf( 'Best Cross Validation Accuracy = %g%% Best c = %g Best g = %g',bestacc,bestc,bestg);
disp(str);
cmd = ['-c ',num2str(bestc),' -g ',num2str(bestg)];
model = svmtrain(train_wine_labels, train_wine, cmd);


%%  仿真测试
p_test= normalized_data3;
t_test = [ones(50,1); 2*ones(160,1); 3*ones(105,1);];
T_sim1 = svmpredict(train_wine_labels,train_wine, model);
T_sim2 = svmpredict(t_test , p_test , model);

%%  数据排序
[train_wine_labels, index_1] = sort(train_wine_labels);
[t_test , index_2] = sort(t_test);

T_sim1 = T_sim1(index_1);
T_sim2 = T_sim2(index_2);
%%
accurate_matches = 0;
total_elements = length(train_wine_labels);
% 遍历每个元素，计算准确匹配率
for i = 1:total_elements
    if train_wine_labels(i) == T_sim1 (i)
        accurate_matches = accurate_matches + 1;
    end
end
% 计算准确率
accuracy1 = accurate_matches / total_elements * 100;
accurate_matches = 0;
total_elements1 = length(t_test);
% 遍历每个元素，计算准确匹配率
for i = 1:total_elements1
    if t_test(i) == T_sim2 (i)
        accurate_matches = accurate_matches + 1;
    end
end
% 计算准确率
accuracy2 = accurate_matches / total_elements1 * 100;



%%  混淆矩阵
figure
cm = confusionchart(train_wine_labels, T_sim1);
cm.Title = '训练集集混淆矩阵';    
figure
cm = confusionchart(t_test, T_sim2);
cm.Title = '测试集混淆矩阵';
