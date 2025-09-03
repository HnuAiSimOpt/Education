

 data1=xlsread("D:\弟弟\只用时域时频域特征\data_x1_PCA.xlsx");
 data2=xlsread("D:\弟弟\只用时域时频域特征\data_x4_PCA.xlsx");
 data3=xlsread("D:\弟弟\只用时域时频域特征\data_x6_PCA.xlsx");
 data11=data1(2:316,2:4);
 data12=data2(2:316,2:4);
 data13=data3(2:316,2:4);
 normalized_data1 = zscore(data11);
 normalized_data2 = zscore(data12);
 normalized_data3 = zscore(data13);
 
P_train = [normalized_data2; normalized_data3; ]';
T_train =[ones(55,1); 2*ones(165,1); 3*ones(95,1);ones(55,1); 2*ones(165,1); 3*ones(95,1)]';
M = size(P_train, 2);

P_test = normalized_data1';
T_test = [ones(40,1); 2*ones(150,1); 3*ones(125,1);]';
N = size(P_test, 2);

%%  数据归一化
[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input );
t_train = T_train;
t_test  = T_test ;

%%  转置以适应模型
train = p_train'; p_test = p_test';
train_label = t_train'; t_test = t_test';

%%  参数设置
cmin=-8;
cmax=8;
gmin=-8;
gmax=8;
v=3;
cstep=1;
gstep=1;
msestep=4.5;
%%  提取最佳参数c和g
[bestCVmse,bestc,bestg] = SVMcgForClass(train_label,train,cmin,cmax,gmin,gmax,v,cstep,gstep,msestep);

%%  建立模型
cmd = [' -c ', num2str(bestc), ' -g ', num2str(bestg)];
model = svmtrain(train_label, train, cmd);


%%  仿真测试
T_sim2 = svmpredict(t_test , p_test , model);
T_sim1 = svmpredict(train_label, train, model);
%%  数据排序
[T_train, index_1] = sort(T_train);
[T_test , index_2] = sort(T_test );

T_sim1 = T_sim1(index_1);
T_sim2 = T_sim2(index_2);
error1 = sum((T_sim1' == T_train)) / M * 100 ;

error2 = sum((T_sim2' == T_test )) / N * 100 ;

%%  混淆矩阵
figure
cm = confusionchart(T_train, T_sim1);
cm.Title = '训练集集混淆矩阵';    
figure
cm = confusionchart(T_test, T_sim2);
cm.Title = '测试集混淆矩阵';
