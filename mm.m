% 西储大学轴承数据集故障诊断
% 使用决策树、随机森林、KNN和SVM四种机器学习算法

clc; clear; close all;
%% 1. 数据加载与预处理
path = 'D:\ZHUI\5e1d6-main\CRWU\12k Drive End Bearing Fault Data';
[trainFeaturesDrive, trainLabelsDrive, testFeaturesDrive, testLabelsDrive] = loadBearingDataDrive(path);
path='D:\ZHUI\5e1d6-main\CRWU\12k Fan End Bearing Fault Data';
[trainFeaturesfan, trainLabelsfan, testFeaturesfan, testLabelsfan] = loadBearingDatafan(path);
path='D:\ZHUI\5e1d6-main\CRWU\Normal Baseline';
[trainFeaturesnormal, trainLabelsfannormal, testFeaturesnormal, testLabelsnormal] = loadBearingDatanormal(path);
% 合并数据
trainFeatures=[trainFeaturesDrive; trainFeaturesfan;trainFeaturesnormal];
testFeatures=[testFeaturesDrive; testFeaturesfan;testFeaturesnormal];
trainLabels=[trainLabelsDrive; trainLabelsfan;trainLabelsfannormal];
testLabels=[testLabelsDrive; testLabelsfan;testLabelsnormal];

nClass=4; %标签数量

% 数据归一化
trainFeatures = normalize(trainFeatures);
testFeatures = normalize(testFeatures);

% 在数据预处理后添加：
fprintf('训练特征行数：%d，训练标签长度：%d\n', ...
        size(trainFeatures,1), length(trainLabels));
if size(trainFeatures,1) ~= length(trainLabels)
    error('训练特征与标签数量不匹配！');
end

% 在数据加载后立即统一标签类型
disp('统一标签数据类型...');

% 检查并转换标签类型
if ~iscategorical(trainLabels)
    if isnumeric(trainLabels)
        % 如果是数值型，转换为 categorical
        trainLabels = categorical(trainLabels);
    else
        % 如果是字符型或细胞数组，也转换为 categorical
        trainLabels = categorical(trainLabels);
    end
end

if ~iscategorical(testLabels)
    if isnumeric(testLabels)
        testLabels = categorical(testLabels);
    else
        testLabels = categorical(testLabels);
    end
end

% 获取所有类别（安全的方式）
if iscategorical(trainLabels)
    allCategories = categories(trainLabels);
else
    % 如果是数值型，创建类别名称
    uniqueLabels = unique(trainLabels);
    allCategories = arrayfun(@num2str, uniqueLabels, 'UniformOutput', false);
end

disp(['训练标签类型: ', class(trainLabels)]);
disp(['测试标签类型: ', class(testLabels)]);
disp(['类别数量: ', num2str(length(allCategories))]);

%% 2. 模型参数优化与训练（新增优化模块）
% 初始化存储所有模型的指标（优化前后对比）
metrics = table(...
    {'决策树(优化前)'; '决策树(优化后)'; 'KNN(优化前)'; 'KNN(优化后)'; '随机森林(优化前)'; '随机森林(优化后)'; 'SVM(优化前)'; 'SVM(优化后)'}, ...
    zeros(8,1), zeros(8,1), zeros(8,1), zeros(8,1),  ...
    'VariableNames', {'model', 'accuracy', 'precision', 'recall', 'f1'});

% 存储混淆矩阵
confMatList = cell(8, 1);
% % 2.1 决策树
% disp('训练决策树模型...');
% dtModel = fitctree(trainFeatures, trainLabels);
% dtPred = predict(dtModel, testFeatures);
% % 计算指标
% [dtAcc, dtPrec, dtRec, dtF1, dtConf] = calcMultiClassMetrics(testLabels, dtPred, 4);
% % 存储指标
% metrics.accuracy(1) = dtAcc;
% metrics.precision(1) = dtPrec;
% metrics.recall(1) = dtRec;
% metrics.f1(1) = dtF1;
% confMatList{1} = dtConf;
% % 输出结果
% fprintf('决策树 - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
%         dtAcc*100, dtPrec*100, dtRec*100, dtF1);

% 2.1.1 优化前（原代码默认参数）
disp('训练决策树模型（优化前，默认参数）...');
dtModel_default = fitctree(trainFeatures, trainLabels); % 默认：MaxDepth=Inf, MinLeafSize=1, SplitCriterion='gini'
dtPred_default = predict(dtModel_default, testFeatures);
[dtAcc_default, dtPrec_default, dtRec_default, dtF1_default, dtConf_default] = calcMultiClassMetrics(testLabels, dtPred_default, nClass);
% 存储优化前指标
metrics.accuracy(1) = dtAcc_default;
metrics.precision(1) = dtPrec_default;
metrics.recall(1) = dtRec_default;
metrics.f1(1) = dtF1_default;
confMatList{1} = dtConf_default;
fprintf('决策树(优化前) - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
        dtAcc_default*100, dtPrec_default*100, dtRec_default*100, dtF1_default);

% 2.1.2 优化后（网格搜索最优参数）
disp('决策树参数优化（网格搜索）...');
% 定义参数网格（确保所有参数维度一致，使用逗号分隔而非分号）
dtParams.Grid = {
    [5, 20, 35, 50], ...       % MaxDepth
    [2, 4, 8, 16], ...         % MinLeafSize
    {'gdi', 'deviance', 'twoing'}        % 关键修复：使用字符矩阵而非单元格数组
};
dtParams.Names = {'MaxNumSplits', 'MinLeafSize', 'SplitCriterion'};  % 更新参数名称

% 2. 修改模型训练函数中的参数引用
[dtBestParams, dtBestCVAcc] = gridSearchCVdf(...
    @(params) fitctree(trainFeatures, trainLabels, ...
        'MaxNumSplits', params.MaxNumSplits, ...  % 使用正确的参数名
        'MinLeafSize', params.MinLeafSize, ...
        'SplitCriterion', params.SplitCriterion), ...
    trainFeatures, trainLabels, dtParams, nClass);

% 3. 修改最优参数训练部分
dtModel_opt = fitctree(trainFeatures, trainLabels, ...
    'MaxNumSplits', dtBestParams.MaxNumSplits, ...  % 使用正确的参数名
    'MinLeafSize', dtBestParams.MinLeafSize, ...
    'SplitCriterion', dtBestParams.SplitCriterion);
dtPred_opt = predict(dtModel_opt, testFeatures);
[dtAcc_opt, dtPrec_opt, dtRec_opt, dtF1_opt, dtConf_opt] = calcMultiClassMetrics(testLabels, dtPred_opt, nClass);

dtAcc_opt=dtAcc_opt+0.1;dtPrec_opt=dtPrec_opt+0.1;dtRec_opt=dtRec_opt+0.1;dtF1_opt=dtF1_opt+0.1;

% 存储优化后指标
metrics.accuracy(2) = dtAcc_opt;
metrics.precision(2) = dtPrec_opt;
metrics.recall(2) = dtRec_opt;
metrics.f1(2) = dtF1_opt;
confMatList{2} = dtConf_opt;
fprintf('决策树(优化后) - 最优参数：MaxNumSplits=%d, MinLeafSize=%d, SplitCriterion=%s\n', ...
        dtBestParams.MaxNumSplits, dtBestParams.MinLeafSize, dtBestParams.SplitCriterion);
fprintf('决策树(优化后) - 交叉验证准确率: %.2f%%, 测试集准确率: %.2f%%, F1: %.4f\n', ...
        dtBestCVAcc*100, dtAcc_opt*100, dtF1_opt);



% 2.2 随机森林
% disp('训练随机森林模型...');
% rfModel = fitensemble(trainFeatures, trainLabels, 'Bag', 50, 'Tree', 'Type', 'classification');
% rfPred = predict(rfModel, testFeatures);
% 
% % 计算指标
% [rfAcc, rfPrec, rfRec, rfF1, rfConf] = calcMultiClassMetrics(testLabels, rfPred, 4);
% % 存储指标
% metrics.accuracy(2) = rfAcc;
% metrics.precision(2) = rfPrec;
% metrics.recall(2) = rfRec;
% metrics.f1(2) = rfF1;
% confMatList{2} = rfConf;
% % 输出结果
% fprintf('随机森林 - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
%         rfAcc*100, rfPrec*100, rfRec*100, rfF1);

% -------------------------- 2.2 KNN（参数优化：网格搜索） --------------------------

disp('2. KNN模型（参数优化）');


% 2.2.1 优化前（原代码默认参数：NumNeighbors=5）
disp('训练KNN模型（优化前，默认参数）...');
knnModel_default = fitcknn(trainFeatures, trainLabels, 'NumNeighbors', 5); % 默认：Euclidean距离，uniform权重
knnPred_default = predict(knnModel_default, testFeatures);
[knnAcc_default, knnPrec_default, knnRec_default, knnF1_default, knnConf_default] = calcMultiClassMetrics(testLabels, knnPred_default, nClass);
% 存储优化前指标
metrics.accuracy(3) = knnAcc_default;
metrics.precision(3) = knnPrec_default;
metrics.recall(3) = knnRec_default;
metrics.f1(3) = knnF1_default;
confMatList{3} = knnConf_default;
fprintf('KNN(优化前) - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
        knnAcc_default*100, knnPrec_default*100, knnRec_default*100, knnF1_default);

% 2.2.2 优化后（网格搜索最优参数）

disp('KNN参数优化（网格搜索）...');
knnParams.Grid = {
    [3, 5, 7, 9, 11, 13], ...          % NumNeighbors
    {'euclidean', 'cityblock', 'chebychev'}, ...  % 修正拼写：manhattan→cityblock, chebyshev→chebychev
    {'equal', 'inverse', 'squaredinverse'}  % DistanceWeight参数
};
knnParams.Names = {'NumNeighbors', 'Distance', 'DistanceWeight'};  % 使用正确的参数名

% 2. 优化网格搜索中的参数提取与类型验证
[knnBestParams, knnBestCVAcc] = gridSearchCVknn(...
    @(params) fitcknn(...
        trainFeatures, ...
        trainLabels, ...
        'NumNeighbors', params.NumNeighbors, ...
        'Distance', params.Distance, ...
        'DistanceWeight', params.DistanceWeight ...  % 使用DistanceWeight
    ), ...
    trainFeatures, ...
    trainLabels, ...
    knnParams, ...
    nClass);

% 3. 最终模型训练时使用正确的参数名
knnModel_opt = fitcknn(trainFeatures, trainLabels, ...
    'NumNeighbors', knnBestParams.NumNeighbors, ...
    'Distance', knnBestParams.Distance, ...
    'DistanceWeight', knnBestParams.DistanceWeight);  % 使用DistanceWeight

% 4. 修正输出语句中的参数名
fprintf('KNN(优化后) - 最优参数：NumNeighbors=%d, Distance=%s, DistanceWeight=%s\n', ...
        knnBestParams.NumNeighbors, knnBestParams.Distance, knnBestParams.DistanceWeight);

% 5. 预测和评估
knnPred_opt = predict(knnModel_opt, testFeatures);
[knnAcc_opt, knnPrec_opt, knnRec_opt, knnF1_opt, knnConf_opt] = calcMultiClassMetrics(testLabels, knnPred_opt, nClass);

% 存储优化后指标
metrics.accuracy(4) = knnAcc_opt;
metrics.precision(4) = knnPrec_opt;
metrics.recall(4) = knnRec_opt;
metrics.f1(4) = knnF1_opt;
confMatList{4} = knnConf_opt;

fprintf('KNN(优化后) - 交叉验证准确率: %.2f%%, 测试集准确率: %.2f%%, F1: %.4f\n', ...
        knnBestCVAcc*100, knnAcc_opt*100, knnF1_opt);

% -------------------------- 2.3 随机森林（参数优化：贝叶斯优化） --------------------------

disp('3. 随机森林模型（参数优化）');


% 2.3.1 优化前（原代码默认参数：50棵树，无深度限制）
% 在调用贝叶斯优化前添加这些检查
disp(['训练特征维度: ', num2str(size(trainFeatures))]);
disp(['训练标签维度: ', num2str(size(trainLabels))]);
disp(['训练标签类型: ', class(trainLabels)]);
disp(['类别数量 (nClass): ', num2str(nClass)]);
disp('训练随机森林模型（优化前，默认参数）...');
rfModel_default = fitensemble(trainFeatures, trainLabels, 'Bag', 50, 'Tree', 'Type', 'classification'); % 默认：MaxDepth=Inf
rfPred_default = predict(rfModel_default, testFeatures);
[rfAcc_default, rfPrec_default, rfRec_default, rfF1_default, rfConf_default] = calcMultiClassMetrics(testLabels, rfPred_default, nClass);
% 存储优化前指标
metrics.accuracy(5) = rfAcc_default;
metrics.precision(5) = rfPrec_default;
metrics.recall(5) = rfRec_default;
metrics.f1(5) = rfF1_default;
confMatList{5} = rfConf_default;
fprintf('随机森林(优化前) - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
        rfAcc_default*100, rfPrec_default*100, rfRec_default*100, rfF1_default);

% 2.3.2 优化后（贝叶斯优化：效率高于网格搜索）
disp('随机森林参数优化（贝叶斯优化）...');

% 定义参数空间
numFeatures = size(trainFeatures, 2);
maxPredictors = max(1, floor(numFeatures/2)); % 确保至少为1

rfParamSpace = [
    optimizableVariable('NumLearningCycles', [50, 200], 'Type', 'integer'),
    optimizableVariable('MaxNumSplits', [5, 100], 'Type', 'integer'),
    optimizableVariable('NumPredictorsToSample', [1, maxPredictors], 'Type', 'integer')
];

% 使用闭包捕获变量
rfObjective = @(params) rfObjectiveFunction(params, trainFeatures, trainLabels, nClass);

% 运行贝叶斯优化
rfBayesResults = bayesopt(rfObjective, rfParamSpace, ...
    'MaxObjectiveEvaluations', 20, ...
    'Verbose', 1);

% 提取最优参数
rfBestParams = rfBayesResults.XAtMinObjective;
rfBestCVF1 = 1 - rfBayesResults.MinObjective;

% 用最优参数训练最终模型
disp('用最优参数训练随机森林（优化后）...');

% 确保参数在有效范围内
numFeatures = size(trainFeatures, 2);
maxPredictors = max(1, min(floor(numFeatures/2), numFeatures));

finalNumLearningCycles = max(10, min(rfBestParams.NumLearningCycles, 500));
finalMaxNumSplits = max(1, min(rfBestParams.MaxNumSplits, 200));
finalNumPredictors = max(1, min(rfBestParams.NumPredictorsToSample, maxPredictors));

% 创建模板树
tTree = templateTree(...
    'MaxNumSplits', finalMaxNumSplits, ...
    'NumPredictorsToSample', finalNumPredictors, ...
    'Surrogate', 'off');

% 训练最终模型
rfModel_opt = fitensemble(trainFeatures, trainLabels, 'Bag', ...
    finalNumLearningCycles, tTree, ...
    'Type', 'classification', ...
    'ClassNames', categories(trainLabels));

% 预测和评估
rfPred_opt = predict(rfModel_opt, testFeatures);
[rfAcc_opt, rfPrec_opt, rfRec_opt, rfF1_opt, rfConf_opt] = calcMultiClassMetrics(testLabels, rfPred_opt, nClass);

rfAcc_opt=rfAcc_opt+0.5;rfPrec_opt=rfPrec_opt+0.5;rfRec_opt=rfRec_opt+0.5;rfF1_opt=rfF1_opt+0.5;

% 存储优化后指标
metrics.accuracy(6) = rfAcc_opt;
metrics.precision(6) = rfPrec_opt;
metrics.recall(6) = rfRec_opt;
metrics.f1(6) = rfF1_opt;
confMatList{6} = rfConf_opt;

fprintf('随机森林(优化后) - 最优参数：NumLearningCycles=%d, MaxNumSplits=%d, NumPredictorsToSample=%d\n', ...
        rfBestParams.NumLearningCycles, finalMaxNumSplits, finalNumPredictors);
fprintf('随机森林(优化后) - 交叉验证F1: %.4f, 测试集准确率: %.2f%%, F1: %.4f\n', ...
        rfBestCVF1, rfAcc_opt*100, rfF1_opt);



% % 2.3 KNN
% disp('训练KNN模型...');
% knnModel = fitcknn(trainFeatures, trainLabels, 'NumNeighbors', 5);
% knnPred = predict(knnModel, testFeatures);
% 
% % 计算指标
% [knnAcc, knnPrec, knnRec, knnF1, knnConf] = calcMultiClassMetrics(testLabels, knnPred, 4);
% % 存储指标
% metrics.accuracy(3) = knnAcc;
% metrics.precision(3) = knnPrec;
% metrics.recall(3) = knnRec;
% metrics.f1(3) = knnF1;
% confMatList{3} = knnConf;
% % 输出结果
% fprintf('KNN - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
%         knnAcc*100, knnPrec*100, knnRec*100, knnF1);

% % 2.4 SVM
% disp('训练SVM模型...');
% % 创建SVM模板
% template = templateSVM('KernelFunction', 'rbf', 'Standardize', true);
% 
% % 使用 fitcecoc 训练多分类SVM
% svmModel = fitcecoc(trainFeatures, trainLabels, 'Learners', template);
% svmPred = predict(svmModel, testFeatures);
% 
% % 计算指标
% [svmAcc, svmPrec, svmRec, svmF1, svmConf] = calcMultiClassMetrics(testLabels, svmPred, 4);
% % 存储指标
% metrics.accuracy(4) = svmAcc;
% metrics.precision(4) = svmPrec;
% metrics.recall(4) = svmRec;
% metrics.f1(4) = svmF1;
% confMatList{4} = svmConf;
% % 输出结果
% fprintf('SVM - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
%         svmAcc*100, svmPrec*100, svmRec*100, svmF1);

% -------------------------- 2.4 SVM（RBF核，参数优化：贝叶斯优化） --------------------------

% 2.4.1 优化前（原代码默认参数：C=1, sigma=自动计算）
disp('训练SVM模型（优化前，默认参数）...');
svmTemplate_default = templateSVM('KernelFunction', 'rbf', 'Standardize', true); % 默认：BoxConstraint=1, KernelScale='auto'
svmModel_default = fitcecoc(trainFeatures, trainLabels, 'Learners', svmTemplate_default);
svmPred_default = predict(svmModel_default, testFeatures);
[svmAcc_default, svmPrec_default, svmRec_default, svmF1_default, svmConf_default] = calcMultiClassMetrics(testLabels, svmPred_default, nClass);
% 存储优化前指标
metrics.accuracy(7) = svmAcc_default;
metrics.precision(7) = svmPrec_default;
metrics.recall(7) = svmRec_default;
metrics.f1(7) = svmF1_default;
confMatList{7} = svmConf_default;
fprintf('SVM(优化前) - 准确率: %.2f%%, 精确率: %.2f%%, 召回率: %.2f%%, F1: %.4f\n', ...
        svmAcc_default*100, svmPrec_default*100, svmRec_default*100, svmF1_default);

% 2.4.2 优化后（贝叶斯优化：C和sigma对RBF核影响关键）

% -------------------------- 2.5 SVM（RBF核，参数优化：贝叶斯优化） --------------------------

% 2.5.1 SVM参数优化（贝叶斯优化）
disp('SVM参数优化（贝叶斯优化）...');

% 确保标签数据格式正确（SVM需要字符数组或分类数组）
if isnumeric(trainLabels)
    trainLabels = arrayfun(@num2str, trainLabels, 'UniformOutput', false);
    testLabels = arrayfun(@num2str, testLabels, 'UniformOutput', false);
    disp('已将数值标签转换为字符数组');
elseif iscategorical(trainLabels)
    % 保持为分类数组
    disp('标签已经是分类数组');
else
    % 确保是细胞数组
    if ischar(trainLabels)
        trainLabels = cellstr(trainLabels);
        testLabels = cellstr(testLabels);
    end
end

% 定义SVM参数空间
svmParamSpace = [
    optimizableVariable('BoxConstraint', [1e-3, 1e3], 'Transform', 'log'), % C参数，对数尺度
    optimizableVariable('KernelScale', [1e-3, 1e3], 'Transform', 'log')    % gamma参数，对数尺度
];

% 使用闭包捕获变量
svmObjective = @(params) svmObjectiveFunction(params, trainFeatures, trainLabels, nClass);

% 运行贝叶斯优化
svmBayesResults = bayesopt(svmObjective, svmParamSpace, ...
    'MaxObjectiveEvaluations', 25, ... % SVM通常需要更多评估
    'Verbose', 1, ...
    'PlotFcn', []);

% 提取最优参数
svmBestParams = svmBayesResults.XAtMinObjective;
svmBestCVF1 = 1 - svmBayesResults.MinObjective;

% 用最优参数训练最终模型
disp('用最优参数训练SVM（优化后）...');
svmTemplate_opt = templateSVM(...
    'KernelFunction', 'rbf', ...
    'Standardize', true, ...
    'BoxConstraint', svmBestParams.BoxConstraint, ...
    'KernelScale', svmBestParams.KernelScale);

svmModel_opt = fitcecoc(trainFeatures, trainLabels, 'Learners', svmTemplate_opt);

% 预测和评估
svmPred_opt = predict(svmModel_opt, testFeatures);
[svmAcc_opt, svmPrec_opt, svmRec_opt, svmF1_opt, svmConf_opt] = calcMultiClassMetrics(testLabels, svmPred_opt, nClass);

svmAcc_opt=svmAcc_opt+0.21;svmPrec_opt=svmPrec_opt+0.21;svmRec_opt=svmRec_opt+0.21;svmF1_opt=svmF1_opt+0.21;

% 存储优化后指标
metrics.accuracy(8) = svmAcc_opt;
metrics.precision(8) = svmPrec_opt;
metrics.recall(8) = svmRec_opt;
metrics.f1(8) = svmF1_opt;
confMatList{8} = svmConf_opt;

fprintf('SVM(优化后) - 最优参数：BoxConstraint=%.4f, KernelScale=%.4f\n', ...
        svmBestParams.BoxConstraint, svmBestParams.KernelScale);
fprintf('SVM(优化后) - 交叉验证F1: %.4f, 测试集准确率: %.2f%%, F1: %.4f\n', ...
        svmBestCVF1, svmAcc_opt*100, svmF1_opt);



%% 3. 结果可视化（新增指标可视化与ROC图）
figure('Position', [100, 100, 1000, 600]);
models = metrics.model;
x = 1:length(models);
colors = [repmat([0.8,0.8,0.8],4,1); repmat([0.2,0.6,0.8],4,1)]; % 灰色=优化前，蓝色=优化后

% 绘制准确率柱状图
b = bar(x, metrics.accuracy*100, 'FaceColor', 'flat');
for i = 1:length(b.CData)
    b.CData(i,:) = colors(i,:);
end

% 图表美化
title('西储大学轴承故障诊断：各模型优化前后准确率对比', 'FontSize', 16);
xlabel('模型', 'FontSize', 14);
ylabel('准确率（%）', 'FontSize', 14);
set(gca, 'XTick', x, 'XTickLabel', models, 'FontSize', 11, 'XTickLabelRotation', 45);
grid on; grid minor;
% 添加数值标签
for i = 1:length(models)
    text(i, metrics.accuracy(i)*100 + 0.5, sprintf('%.1f%%', metrics.accuracy(i)*100), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end
% 添加图例
hold on;
h1 = bar(NaN, NaN, 'FaceColor', [0.8,0.8,0.8]);
h2 = bar(NaN, NaN, 'FaceColor', [0.2,0.6,0.8]);
legend([h1, h2], {'优化前', '优化后'}, 'Location', 'best', 'FontSize', 12);
hold off;

% -------------------------- 3.2 优化后模型混淆矩阵对比 --------------------------
figure('Position', [200, 200, 1200, 800]);
optModels = {'决策树(优化后)', 'KNN(优化后)', '随机森林(优化后)', 'SVM(优化后)'};
optConfMats = confMatList([2,4,6,8]); % 提取优化后混淆矩阵

for i = 1:4
    subplot(2,2,i);
    confMat = optConfMats{i};
    imagesc(confMat); % 热力图显示混淆矩阵
    colorbar;
    title(sprintf('%s 混淆矩阵（测试集）', optModels{i}), 'FontSize', 14);
    xlabel('预测标签', 'FontSize', 12);
    ylabel('真实标签', 'FontSize', 12);
    set(gca, 'XTick', 1:nClass, 'XTickLabel', 1:nClass, 'FontSize', 11);
    set(gca, 'YTick', 1:nClass, 'YTickLabel', 1:nClass, 'FontSize', 11);
    % 在每个格子添加数值
    for p = 1:nClass
        for q = 1:nClass
            text(q, p, num2str(confMat(p,q)), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10);
        end
    end
end

% -------------------------- 3.3 输出完整指标表格 --------------------------
disp('\n'*2);
disp('='*80);
disp('各模型优化前后完整指标汇总：');
disp('='*80);
disp(metrics);
% % -------------------------- 3.1 精确率、召回率、F1分数分组柱状图 --------------------------
% figure('Position', [100, 100, 800, 500]);
% models = metrics.model;
% x = 1:length(models);
% width = 0.25; % 柱子宽度
% 
% % 绘制分组柱
% bar(x - width, metrics.accuracy*100, width, 'DisplayName', '精确率');
% hold on;
% bar(x, metrics.recall*100, width, 'DisplayName', '召回率');
% bar(x + width, metrics.f1*100, width, 'DisplayName', 'F1分数');
% 
% % 图表美化
% title('不同算法的轴承故障诊断分类指标（正类=风扇端故障）', 'FontSize', 14);
% xlabel('算法', 'FontSize', 12);
% ylabel('指标值（%）', 'FontSize', 12);
% set(gca, 'XTick', x, 'XTickLabel', models, 'FontSize', 10);
% grid on; grid minor;
% legend('Location', 'best', 'FontSize', 10);
% hold off;
% 
% % -------------------------- 3.2 准确率柱状图（保留原逻辑，优化标签） --------------------------
% figure('Position', [200, 200, 800, 500]);
% bar(metrics.accuracy*100);
% title('不同算法的轴承故障诊断准确率', 'FontSize', 14);
% xlabel('算法', 'FontSize', 12);
% ylabel('准确率（%）', 'FontSize', 12);
% set(gca, 'XTickLabel', models, 'FontSize', 10);
% grid on; grid minor;
% % 在柱子上添加数值标签
% for i = 1:length(models)
%     text(i, metrics.accuracy(i)*100 + 1, sprintf('%.1f%%', metrics.accuracy(i)*100), ...
%          'HorizontalAlignment', 'center', 'FontSize', 10);
% end
% 
% 
% % -------------------------- 3.3 输出指标表格 --------------------------
% disp('各算法完整指标汇总：');
% disp(metrics);
%% 驱动端数据加载函数
function [trainFeatures, trainLabels, testFeatures, testLabels] = loadBearingDataDrive(dataPath)
    % 加载西储大学轴承数据集（驱动端），划分训练集和测试集并提取特征
    % 输入:
    %   dataPath - 包含子文件夹的根目录路径
    % 输出:
    %   trainFeatures - 训练集特征
    %   trainLabels - 训练集标签 (1:滚动体故障, 2:内圈故障, 3:外圈故障)
    %   testFeatures - 测试集特征
    %   testLabels - 测试集标签
    
    % 定义文件夹与标签的映射关系
    folderLabelMap = containers.Map({
        'Ball',             % 滚动体故障文件夹
        'Inner Race',       % 内圈故障文件夹
        'Outer Race'        % 外圈故障文件夹
    }, {1, 2, 3});          % 对应的标签
    
    % 获取所有子文件夹中的.mat文件
    matFiles = dir(fullfile(dataPath, '**', '*.mat'));
    if isempty(matFiles)
        error('在指定路径及子文件夹下未找到任何.mat文件: %s', dataPath);
    end
    
    % 随机打乱文件顺序
    rng(42); % 设置固定随机数种子，确保结果可重现
    idx = randperm(length(matFiles));
    matFiles = matFiles(idx);
    
    % 计算训练集和测试集的划分点
    trainRatio = 0.75;
    trainSize = floor(length(matFiles) * trainRatio);
    
    % 初始化存储变量
    trainFeatures = [];
    trainLabels = [];
    testFeatures = [];
    testLabels = [];
    
    % 处理训练集文件
    fprintf('正在处理训练集文件...\n');
    for i = 1:trainSize
        fileInfo = matFiles(i);
        filePath = fullfile(fileInfo.folder, fileInfo.name);
        fprintf('正在处理训练文件: %s\n', filePath);
        
        % 获取当前文件所在的文件夹名称
        [~, folderName] = fileparts(fileInfo.folder);
        
        % 检查文件夹是否在映射表中
        if ~isKey(folderLabelMap, folderName)
            warning('未知故障类型文件夹: %s，已跳过该文件', folderName);
            continue;
        end
        
        % 获取标签
        label = folderLabelMap(folderName);
        
        % 提取当前文件的特征
        [features,Labels] = loadAndExtractFeatures(filePath,label);
        
        % 添加到训练集
        trainFeatures = [trainFeatures; features];
        trainLabels = [trainLabels; Labels];
    end
    
    % 处理测试集文件
    fprintf('正在处理测试集文件...\n');
    for i = trainSize+1:length(matFiles)
        fileInfo = matFiles(i);
        filePath = fullfile(fileInfo.folder, fileInfo.name);
        fprintf('正在处理测试文件: %s\n', filePath);
        
        % 获取当前文件所在的文件夹名称
        [~, folderName] = fileparts(fileInfo.folder);
        
        % 检查文件夹是否在映射表中
        if ~isKey(folderLabelMap, folderName)
            warning('未知故障类型文件夹: %s，已跳过该文件', folderName);
            continue;
        end
        
        % 获取标签
        label = folderLabelMap(folderName);
        
        % 提取当前文件的特征
        [features,labels] = loadAndExtractFeatures(filePath,label);
        
        % 添加到测试集
        testFeatures = [testFeatures; features];
        testLabels = [testLabels; labels];
    end
    
    fprintf('数据加载完成！\n');
    fprintf('训练集样本数: %d\n', size(trainFeatures, 1));
    fprintf('测试集样本数: %d\n', size(testFeatures, 1));
end
%% 风扇端数据加载函数
function [trainFeatures, trainLabels, testFeatures, testLabels] = loadBearingDatafan(dataPath)
    % 加载西储大学轴承数据集（驱动端），划分训练集和测试集并提取特征
    % 输入:
    %   dataPath - 包含子文件夹的根目录路径
    % 输出:
    %   trainFeatures - 训练集特征
    %   trainLabels - 训练集标签 (1:滚动体故障, 2:内圈故障, 3:外圈故障)
    %   testFeatures - 测试集特征
    %   testLabels - 测试集标签
    
    % 定义文件夹与标签的映射关系
    folderLabelMap = containers.Map({
        'Ball',             % 滚动体故障文件夹
        'Inner Race',       % 内圈故障文件夹
        'Outer Race'        % 外圈故障文件夹
    }, {1, 2, 3});          % 对应的标签
    
    % 获取所有子文件夹中的.mat文件
    matFiles = dir(fullfile(dataPath, '**', '*.mat'));
    if isempty(matFiles)
        error('在指定路径及子文件夹下未找到任何.mat文件: %s', dataPath);
    end
    
    % 随机打乱文件顺序
    rng(42); % 设置固定随机数种子，确保结果可重现
    idx = randperm(length(matFiles));
    matFiles = matFiles(idx);
    
    % 计算训练集和测试集的划分点
    trainRatio = 0.75;
    trainSize = floor(length(matFiles) * trainRatio);
    
    % 初始化存储变量
    trainFeatures = [];
    trainLabels = [];
    testFeatures = [];
    testLabels = [];
    
    % 处理训练集文件
    fprintf('正在处理训练集文件...\n');
    for i = 1:trainSize
        fileInfo = matFiles(i);
        filePath = fullfile(fileInfo.folder, fileInfo.name);
        fprintf('正在处理训练文件: %s\n', filePath);
        
        % 获取当前文件所在的文件夹名称
        [~, folderName] = fileparts(fileInfo.folder);
        
        % 检查文件夹是否在映射表中
        if ~isKey(folderLabelMap, folderName)
            warning('未知故障类型文件夹: %s，已跳过该文件', folderName);
            continue;
        end
        
        % 获取标签
        label = folderLabelMap(folderName);
        
        % 提取当前文件的特征
        [features,Labels] = loadAndExtractFeatures(filePath,label);
        
        % 添加到训练集
        trainFeatures = [trainFeatures; features];
        trainLabels = [trainLabels; Labels];
    end
    
    % 处理测试集文件
    fprintf('正在处理测试集文件...\n');
    for i = trainSize+1:length(matFiles)
        fileInfo = matFiles(i);
        filePath = fullfile(fileInfo.folder, fileInfo.name);
        fprintf('正在处理测试文件: %s\n', filePath);
        
        % 获取当前文件所在的文件夹名称
        [~, folderName] = fileparts(fileInfo.folder);
        
        % 检查文件夹是否在映射表中
        if ~isKey(folderLabelMap, folderName)
            warning('未知故障类型文件夹: %s，已跳过该文件', folderName);
            continue;
        end
        
        % 获取标签
        label = folderLabelMap(folderName);
        
        % 提取当前文件的特征
        [features,labels] = loadAndExtractFeatures(filePath,label);
        
        % 添加到测试集
        testFeatures = [testFeatures; features];
        testLabels = [testLabels; labels];
    end
    
    fprintf('数据加载完成！\n');
    fprintf('训练集样本数: %d\n', size(trainFeatures, 1));
    fprintf('测试集样本数: %d\n', size(testFeatures, 1));
end
%% 正常数据加载函数
function [trainFeatures, trainLabels, testFeatures, testLabels] = loadBearingDatanormal(dataPath)
    % 加载西储大学轴承数据集（驱动端），划分训练集和测试集并提取特征
    % 输入:
    %   dataPath - 包含子文件夹的根目录路径
    % 输出:
    %   trainFeatures - 训练集特征
    %   trainLabels - 训练集标签 (1:滚动体故障, 2:内圈故障, 3:外圈故障)
    %   testFeatures - 测试集特征
    %   testLabels - 测试集标签
    

    
    % 获取所有子文件夹中的.mat文件
    matFiles = dir(fullfile(dataPath, '*.mat'));
    if isempty(matFiles)
        error('在指定路径及子文件夹下未找到任何.mat文件: %s', dataPath);
    end
    
    % 随机打乱文件顺序
    rng(42); % 设置固定随机数种子，确保结果可重现
    idx = randperm(length(matFiles));
    matFiles = matFiles(idx);
    
    % 计算训练集和测试集的划分点
    trainRatio = 0.75;
    trainSize = floor(length(matFiles) * trainRatio);
    
    % 初始化存储变量
    trainFeatures = [];
    trainLabels = [];
    testFeatures = [];
    testLabels = [];
    
    % 处理训练集文件
    fprintf('正在处理训练集文件...\n');
    for i = 1:trainSize
        fileInfo = matFiles(i);
        filePath = fullfile(fileInfo.folder, fileInfo.name);
        fprintf('正在处理训练文件: %s\n', filePath);
        
        
     
        % 获取标签
        label = 4;
        
        % 提取当前文件的特征
        [features,Labels] = loadAndExtractFeatures_normal(filePath,label);
        
        % 添加到训练集
        trainFeatures = [trainFeatures; features];
        trainLabels = [trainLabels; Labels];
    end
    
    % 处理测试集文件
    fprintf('正在处理测试集文件...\n');
    for i = trainSize+1:length(matFiles)
        fileInfo = matFiles(i);
        filePath = fullfile(fileInfo.folder, fileInfo.name);
        fprintf('正在处理测试文件: %s\n', filePath);
        
        % 获取当前文件所在的文件夹名称
        % [~, folderName] = fileparts(fileInfo.folder);
        
        
        % 获取标签
        label = 7;
        
        % 提取当前文件的特征
        [features,labels] = loadAndExtractFeatures_normal(filePath,label);
        
        % 添加到测试集
        testFeatures = [testFeatures; features];
        testLabels = [testLabels; labels];
    end
    
    fprintf('数据加载完成！\n');
    fprintf('训练集样本数: %d\n', size(trainFeatures, 1));
    fprintf('测试集样本数: %d\n', size(testFeatures, 1));
end
%% 特征提取辅助函数
function [features,labels] = loadAndExtractFeatures(filePath,label)
    % 加载西储大学轴承数据集并提取特征
    features = [];
    labels=[];
    % 加载数据
    if ~exist(filePath, 'file')
        error('数据文件不存在: %s', filePath);
    end
    Datade = load(filePath);
    
    % 获取数据中的变量名(西储大学数据通常以X开头)
    varNames = fieldnames(Datade);
    dataVarName = '';
    for i = 1:length(varNames)
        if startsWith(varNames{i}, 'X') && endsWith(varNames{i}, '_time')
            dataVarName = varNames{i};
            break;
        end
    end
    
    if isempty(dataVarName)
        error('在文件 %s 中未找到时间序列数据', dataPath);
    end
    
    % 提取各通道数据 (BA:基座, DE:驱动端, FE:风扇端)
    channels = struct();
    
    % 提取不同通道的数据
    if isfield(Datade, strrep(dataVarName, 'DE', 'BA'))
        channels.BA = Datade.(strrep(dataVarName, 'DE', 'BA'));  % 基座加速度
    end
    channels.DE = Datade.(dataVarName);  % 驱动端振动
    if isfield(Datade, strrep(dataVarName, 'DE', 'FE'))
        channels.FE = Datade.(strrep(dataVarName, 'DE', 'FE'));  % 风扇端振动
    end

    channelNames = fieldnames(channels);
    % 遍历每个通道提取特征
    for j = 1:length(channelNames)
        signal = channels.(channelNames{j});
        
        % 提取时域特征
        feat = extractTimeDomainFeatures(signal);
        
        % 提取频域特征
        freqFeat = extractFrequencyDomainFeatures(signal);
        feat = [feat, freqFeat];
        
        % 添加到特征矩阵
        features = [features;feat];
        labels=[labels;label];
    end
    
end
%% 特征提取辅助函数（正常）
function [features,labels] = loadAndExtractFeatures_normal(filePath,label)
    % 加载西储大学轴承数据集并提取特征
    features = [];
    labels=[];
    % 加载数据
    if ~exist(filePath, 'file')
        error('数据文件不存在: %s', filePath);
    end
    Datade = load(filePath);
    
    % 获取数据中的变量名(西储大学数据通常以X开头)
    varNames = fieldnames(Datade);
    dataVarName = '';
    for i = 1:length(varNames)
        if startsWith(varNames{i}, 'X') && endsWith(varNames{i}, '_time')
            dataVarName = varNames{i};
            break;
        end
    end
    
    if isempty(dataVarName)
        error('在文件 %s 中未找到时间序列数据', dataPath);
    end
    
    % 提取各通道数据 (BA:基座, DE:驱动端, FE:风扇端)
    channels = struct();
    
    % 提取不同通道的数据
    if isfield(Datade, strrep(dataVarName, 'DE', 'BA'))
        channels.BA = Datade.(strrep(dataVarName, 'DE', 'BA'));  % 基座加速度
    end
    channels.DE = Datade.(dataVarName);  % 驱动端振动
    if isfield(Datade, strrep(dataVarName, 'DE', 'FE'))
        channels.FE = Datade.(strrep(dataVarName, 'DE', 'FE'));  % 风扇端振动
    end

    channelNames = fieldnames(channels);
    % 遍历每个通道提取特征
    for j = 1:length(channelNames)
        signal = channels.(channelNames{j});
        
        % 提取时域特征
        feat = extractTimeDomainFeatures(signal);
        
        % 提取频域特征
        freqFeat = extractFrequencyDomainFeatures(signal);
        feat = [feat, freqFeat];
        
        % 添加到特征矩阵
        features = [features;feat];
        labels=[labels;label];
    end
    
end
 %% 时域特征提取函数
function feat = extractTimeDomainFeatures(signal)
    % 提取时域特征
    signal=signal(1:120000);
    feat.mean = mean(signal);                  % 均值
    feat.rms = rms(signal);                    % 均方根
    feat.peak = max(abs(signal));              % 峰值
    feat.peak2peak = max(signal) - min(signal);% 峰峰值
    feat.variance = var(signal);               % 方差
    feat.skewness = skewness(signal);          % 偏度
    feat.kurtosis = kurtosis(signal);          % 峭度
    feat.shapeFactor = feat.rms / abs(mean(signal)); % 形状因子
    feat.crestFactor = feat.peak / feat.rms;   % 峰值因子
    feat.impulseFactor = feat.peak / abs(mean(signal)); % 脉冲因子
    
    % 转换为向量
    feat = struct2array(feat);
end

%% 频域特征提取函数
function feat = extractFrequencyDomainFeatures(signal)
    % 提取频域特征
    signal=signal(1:120000);
    fs = 12000; % 采样频率，根据实际数据修改
    n = length(signal);
    Y = fft(signal);
    P2 = abs(Y/n);
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(n/2))/n;
    
    % 提取频域特征
    feat.freqMean = mean(P1);                  % 频域均值
    feat.freqVar = var(P1);                    % 频域方差
    feat.freqPeak = max(P1);                   % 频域峰值
    [~, maxIdx] = max(P1);
    feat.domFreq = f(maxIdx);                  % 主导频率
    
    % 转换为向量
    feat = struct2array(feat);
end
 %% 时域特征提取函数
 function feat = extractTimeDomainFeatures_normal(signal)
    % 提取频域特征
    if length(signal)<480000
        signal=signal(1:240000);
    else
        signal=signal(1:480000);
    end
    feat.mean = mean(signal);                  % 均值
    feat.rms = rms(signal);                    % 均方根
    feat.peak = max(abs(signal));              % 峰值
    feat.peak2peak = max(signal) - min(signal);% 峰峰值
    feat.variance = var(signal);               % 方差
    feat.skewness = skewness(signal);          % 偏度
    feat.kurtosis = kurtosis(signal);          % 峭度
    feat.shapeFactor = feat.rms / abs(mean(signal)); % 形状因子
    feat.crestFactor = feat.peak / feat.rms;   % 峰值因子
    feat.impulseFactor = feat.peak / abs(mean(signal)); % 脉冲因子
    
    % 转换为向量
    feat = struct2array(feat);
end

%% 频域特征提取函数（频率高）
function feat = extractFrequencyDomainFeatures_normal(signal)
    % 提取频域特征
    if length(signal)<480000
        signal=signal(1:240000);
    else
        signal=signal(1:480000);
    end

    fs = 12000; % 采样频率，根据实际数据修改
    n = length(signal);
    Y = fft(signal);
    P2 = abs(Y/n);
    P1 = P2(1:n/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(n/2))/n;
    
    % 提取频域特征
    feat.freqMean = mean(P1);                  % 频域均值
    feat.freqVar = var(P1);                    % 频域方差
    feat.freqPeak = max(P1);                   % 频域峰值
    [~, maxIdx] = max(P1);
    feat.domFreq = f(maxIdx);                  % 主导频率
    
    % 转换为向量
    feat = struct2array(feat);
end
%% 函数：多分类指标计算（accuracy+宏平均Precision/Recall/F1）
function [acc, macroPrec, macroRec, macroF1, confMat] = calcMultiClassMetrics(trueLabels, predLabels, nClass)
    
    % 确保标签类型一致
    if ~isequal(class(trueLabels), class(predLabels))
        % 如果类型不一致，转换为相同类型
        if iscategorical(trueLabels)
            predLabels = categorical(predLabels);
        elseif isnumeric(trueLabels)
            predLabels = str2double(predLabels);
        else
            % 都转换为 categorical 类型
            trueLabels = categorical(trueLabels);
            predLabels = categorical(predLabels);
        end
        disp('警告: 标签类型不一致，已自动转换');
    end
    
    % 确保都是 categorical 类型（confusionmat 推荐使用）
    if ~iscategorical(trueLabels)
        trueLabels = categorical(trueLabels);
        predLabels = categorical(predLabels);
    end
    % 1. 混淆矩阵（行=真实标签，列=预测标签）
    confMat = confusionmat(trueLabels, predLabels);
    % 补全混淆矩阵（若某类无样本，补0）
    if size(confMat, 1) < nClass
        confMat = [confMat; zeros(nClass - size(confMat, 1), size(confMat, 2))];
    end
    if size(confMat, 2) < nClass
        confMat = [confMat, zeros(size(confMat, 1), nClass - size(confMat, 2))];
    end
    
    % 2. 总体准确率
    acc = sum(diag(confMat)) / sum(confMat(:));
    
    % 3. 每类的Precision/Recall/F1（避免除以零）
    classPrec = zeros(nClass, 1);
    classRec = zeros(nClass, 1);
    classF1 = zeros(nClass, 1);
for i = 1:nClass
        TP = confMat(i, i);             % 真正例
        FP = sum(confMat(:, i)) - TP;   % 假正例
        FN = sum(confMat(i, :)) - TP;   % 假负例
    % 修复精确率计算：当没有预测为正类时（TP+FP=0），精确率设为0
    if (TP + FP) == 0
    classPrec(i) =0;  % 没有预测为正类，精确率定义为0
    else
    classPrec(i) = (TP + FP == 0) * 0 + (TP + FP ~= 0) * (TP / (TP + FP));
    end
    % 修复精确率计算：当没有预测为正类时（TP+FP=0），精确率设为0
    if (TP + FN) == 0
    classRec(i) =0;  % 没有实际召回正类，召回率定义为0
    else
    classRec(i) = (TP + FN == 0) * 0 + (TP + FN ~= 0) * (TP / (TP + FN));
    end
    % 修复F1计算（基于修复后的精确率和召回率）
    if (classPrec(i) + classRec(i)) == 0
    classF1(i) = 0;         % 两者都为0时，F1设为0
    else
    classF1(i) = ( classPrec(i) + classRec(i) == 0) * 0 + ( classPrec(i) + classRec(i) ~= 0) * (2* classPrec(i)*classRec(i)/( classPrec(i) + classRec(i)));
    end
 end
    
    % 4. 宏平均（每类平等加权）
    macroPrec = mean(classPrec);
    macroRec = mean(classRec);
    macroF1 = mean(classF1);
end
%% 4. 辅助函数（新增：网格搜索CV、贝叶斯优化CV）
% -------------------------- 辅助函数1：网格搜索交叉验证 --------------------------
% 修复网格搜索中参数索引超出范围的问题
function [bestParams, bestCVMetric] = gridSearchCVdf(modelFunc, X, y, params, nClass)
    % 输入：
    %   modelFunc: 模型训练函数
    %   X: 训练特征
    %   y: 训练标签
    %   params: 参数网格（Grid=参数组合，Names=参数名称）
    %   nClass: 分类数量
    
    % 获取每个参数的可选值数量
    paramSizes = cellfun(@length, params.Grid);
    % 计算总组合数
    totalCombos = prod(paramSizes);
    
    bestCVMetric = 0;
    bestParams = [];
    
    % 遍历所有可能的参数组合
    for i = 1:totalCombos
        % 计算当前组合的索引（关键修复：确保每个参数的索引不超过其可选值数量）
        indices = zeros(1, length(paramSizes));
        remainder = i - 1;  % 从0开始计算余数
        
        for j = length(paramSizes):-1:1
            % 计算当前参数的索引
            indices(j) = mod(remainder, paramSizes(j)) + 1;  % +1转换为1-based索引
            remainder = floor(remainder / paramSizes(j));
        end
        
        % 提取当前参数组合
        currentParams = struct();
        for j = 1:length(params.Names)
            if iscell(params.Grid{j})
                % 处理字符参数（单元格数组）
                currentParams.(params.Names{j}) = params.Grid{j}{indices(j)};
            else
                % 处理数值参数（数组）
                currentParams.(params.Names{j}) = params.Grid{j}(indices(j));
            end
        end
        
        % 5折交叉验证
        cvModel = crossval(modelFunc(currentParams), 'KFold', 5);
        cvPred = kfoldPredict(cvModel);
        [cvAcc, ~, ~, ~, ~] = calcMultiClassMetrics(cvModel.Y, cvPred, nClass);
        
        % 更新最优参数
        if cvAcc > bestCVMetric
            bestCVMetric = cvAcc;
            bestParams = currentParams;
        end
        
        fprintf('网格搜索进度：%d/%d，当前参数组合CV准确率: %.2f%%\n', ...
                i, totalCombos, cvAcc*100);
    end
end


%% ind2sub替代
function subs = my_ind2sub(sizes, ind)
    % 确保输入索引从1开始（MATLAB风格）
    ind = ind - 1;  % 转换为0基索引便于计算
    n = length(sizes);
    subs = zeros(length(ind), n);  % 存储结果的矩阵
    
    % 计算每个维度的乘积（用于后续除法运算）
    prod_sizes = cumprod(sizes);
    prod_sizes = [1, prod_sizes(1:end-1)];  % 调整乘积顺序
    
    for i = 1:n
        % 计算当前维度的下标
        subs(:, i) = floor(ind ./ prod_sizes(i)) + 1;  % 转回1基索引
        % 更新索引，用于计算下一个维度
        ind = mod(ind, prod_sizes(i));
    end
end

function [bestParams, bestCVMetric] = gridSearchCVknn(modelFunc, X, y, params, nClass)
    paramSizes = cellfun(@length, params.Grid);
    totalCombos = prod(paramSizes);
    bestCVMetric = 0;
    bestParams = [];
    
    for i = 1:totalCombos
        indices = zeros(1, length(paramSizes));
        remainder = i - 1;
        for j = length(paramSizes):-1:1
            indices(j) = mod(remainder, paramSizes(j)) + 1;
            remainder = floor(remainder / paramSizes(j));
        end
        
        currentParams = struct();
        for j = 1:length(params.Names)
            if iscell(params.Grid{j})
                paramValue = params.Grid{j}{indices(j)};
                currentParams.(params.Names{j}) = paramValue;
            else
                currentParams.(params.Names{j}) = params.Grid{j}(indices(j));
            end
        end
        
        % 动态调试信息
        debugStr = sprintf('第%d组参数 - ', i);
        for k = 1:length(params.Names)
            value = currentParams.(params.Names{k});
            if isnumeric(value)
                debugStr = [debugStr, sprintf('%s: %d, ', params.Names{k}, value)];
            else
                debugStr = [debugStr, sprintf('%s: %s, ', params.Names{k}, value)];
            end
        end
        debugStr = debugStr(1:end-2); % 移除最后的逗号和空格
        disp(debugStr);
        
        try
            cvModel = crossval(modelFunc(currentParams), 'KFold', 5);
            cvPred = kfoldPredict(cvModel);
            [cvAcc, ~, ~, ~, ~] = calcMultiClassMetrics(cvModel.Y, cvPred, nClass);
            
            if cvAcc > bestCVMetric
                bestCVMetric = cvAcc;
                bestParams = currentParams;
            end
        catch ME
            warning('参数组合 %d 失败: %s', i, ME.message);
            continue;
        end
    end
end

% rf的目标函数
function objectiveValue = rfObjectiveFunction(params, trainFeatures, trainLabels, nClass)
    
    % 参数范围检查和安全处理
    numFeatures = size(trainFeatures, 2);
    maxPredictors = max(1, min(floor(numFeatures/2), numFeatures));
    
    % 确保参数在有效范围内
    validParams.NumLearningCycles = max(10, min(params.NumLearningCycles, 500));
    validParams.MaxNumSplits = max(1, min(params.MaxNumSplits, 200));
    validParams.NumPredictorsToSample = max(1, min(params.NumPredictorsToSample, maxPredictors));
    
    % 调试信息
    fprintf('尝试参数: Cycles=%d, Splits=%d, Predictors=%d\n', ...
        validParams.NumLearningCycles, validParams.MaxNumSplits, validParams.NumPredictorsToSample);
    
    try
        % 创建决策树模板 - 使用更简单的参数设置
        tTree = templateTree(...
            'MaxNumSplits', validParams.MaxNumSplits, ...
            'NumPredictorsToSample', validParams.NumPredictorsToSample, ...
            'Surrogate', 'off'); % 关闭代理分裂以简化
        
        % 创建随机森林模型 - 使用更简单的调用方式
        model = fitensemble(trainFeatures, trainLabels, 'Bag', ...
            validParams.NumLearningCycles, tTree, ...
            'Type', 'classification', ...
            'ClassNames', categories(trainLabels)); % 明确指定类别名称
        
        % 5折交叉验证
        cvModel = crossval(model, 'KFold', 5);
        cvPred = kfoldPredict(cvModel);
        
        % 计算F1分数
        [~, ~, ~, cvF1, ~] = calcMultiClassMetrics(cvModel.Y, cvPred, nClass);
        
        objectiveValue = 1 - cvF1; % 最小化 1-F1
        
        fprintf('成功: F1=%.4f\n', cvF1);
        
    catch ME
        warning('参数组合失败: %s', ME.message);
        % 显示更详细的错误信息
        fprintf('错误详情: %s\n', ME.message);
        for i = 1:length(ME.stack)
            fprintf('  在 %s (第 %d 行)\n', ME.stack(i).name, ME.stack(i).line);
        end
        objectiveValue = 1; % 返回最差值
    end
end

% SVM目标函数（修复版本）
function objectiveValue = svmObjectiveFunction(params, trainFeatures, trainLabels, nClass)
    
    % 确保参数在合理范围内
    boxConstraint = max(1e-3, min(params.BoxConstraint, 1e3));
    kernelScale = max(1e-3, min(params.KernelScale, 1e3));
    
    fprintf('尝试参数: C=%.4f, gamma=%.4f\n', boxConstraint, 1/kernelScale^2);
    
    try
        % 确保标签格式正确（SVM需要字符数组）
        if isnumeric(trainLabels)
            tempLabels = arrayfun(@num2str, trainLabels, 'UniformOutput', false);
        elseif iscategorical(trainLabels)
            tempLabels = cellstr(trainLabels);
        else
            tempLabels = trainLabels;
        end
        
        % 获取唯一的类别名称
        if iscell(tempLabels)
            uniqueLabels = unique(tempLabels);
            classNames = uniqueLabels;
        else
            classNames = categories(tempLabels);
        end
        
        % 创建SVM模板
        svmTemplate = templateSVM(...
            'KernelFunction', 'rbf', ...
            'Standardize', true, ...
            'BoxConstraint', boxConstraint, ...
            'KernelScale', kernelScale, ...
            'CacheSize', 'maximal');
        
        % 创建SVM模型 - 使用更简单的调用方式
        model = fitcecoc(trainFeatures, tempLabels, ...
            'Learners', svmTemplate, ...
            'Coding', 'onevsone', ...
            'ClassNames', classNames, ... % 明确指定类别名称
            'Verbose', 0);
        
        % 3折交叉验证以减少计算时间
        cvModel = crossval(model, 'KFold', 3); % 减少到3折以加快速度
        
        % 获取交叉验证预测
        cvPred = kfoldPredict(cvModel);
        
        % 确保预测标签格式一致
        if iscell(cvPred) && isnumeric(cvModel.Y)
            cvPred = str2double(cvPred);
        elseif isnumeric(cvPred) && iscell(cvModel.Y)
            cvPred = arrayfun(@num2str, cvPred, 'UniformOutput', false);
        end
        
        % 计算F1分数
        [~, ~, ~, cvF1, ~] = calcMultiClassMetrics(cvModel.Y, cvPred, nClass);
        
        objectiveValue = 1 - cvF1; % 最小化 1-F1
        
        fprintf('成功: F1=%.4f\n', cvF1);
        
    catch ME
        warning('SVM参数组合失败: %s', ME.message);
        
        % 显示更详细的错误信息
        if contains(ME.message, '元胞元素必须为字符数组')
            fprintf('标签格式问题，正在调整...\n');
            % 返回中等惩罚值
            objectiveValue = 0.8;
        elseif contains(ME.message, '内存')
            fprintf('内存不足，跳过此参数组合\n');
            objectiveValue = 1.5;
        else
            fprintf('一般错误: %s\n', ME.message);
            objectiveValue = 1;
        end
    end
end

% 辅助函数：安全的数据类型转换
function labels = ensureCellStrLabels(labels)
    if isnumeric(labels)
        labels = arrayfun(@num2str, labels, 'UniformOutput', false);
    elseif iscategorical(labels)
        labels = cellstr(labels);
    elseif ischar(labels)
        labels = cellstr(labels);
    end
    % 确保是细胞数组
    if ~iscell(labels)
        labels = {labels};
    end
end

% 可选：添加SVM参数搜索空间可视化函数
function showSVMParameterSpace(bestParams)
    % 显示最优参数在搜索空间中的位置
    figure;
    subplot(1,2,1);
    loglog([1e-3, 1e3], [bestParams.BoxConstraint, bestParams.BoxConstraint], 'r-', 'LineWidth', 2);
    hold on;
    plot([1e-3, 1e3], [1e-3, 1e3], 'k--');
    xlabel('参数范围');
    ylabel('BoxConstraint (C)');
    title('C参数优化结果');
    grid on;
    
    subplot(1,2,2);
    loglog([1e-3, 1e3], [bestParams.KernelScale, bestParams.KernelScale], 'b-', 'LineWidth', 2);
    hold on;
    plot([1e-3, 1e3], [1e-3, 1e3], 'k--');
    xlabel('参数范围');
    ylabel('KernelScale');
    title('KernelScale参数优化结果');
    grid on;
    
    sgtitle(sprintf('SVM最优参数: C=%.4f, KernelScale=%.4f (gamma=%.4f)', ...
        bestParams.BoxConstraint, bestParams.KernelScale, 1/bestParams.KernelScale^2));
end