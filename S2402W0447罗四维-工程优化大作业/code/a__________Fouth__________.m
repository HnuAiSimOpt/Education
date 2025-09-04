clc;
clear;
close all;
addpath(genpath('PR/'))
%-------------代理模型用于多目标优化算例-------------%
theta = [0.01 0.01];
lob = [1e-2 1e-2];
upb = [0.1 0.1];
variables_lob = [3 28];
variables_upb = [4 33];
range = variables_upb - variables_lob;
Decision_variables = importdata('x_1.txt');
testX = importdata('testX.txt');
testY1 = importdata('testY1.txt');
testY2 = importdata('testY2.txt');
Decision_variables =...
    (Decision_variables - variables_lob) ./ repmat(range, size(Decision_variables, 1), 1);

testX =(testX - variables_lob) ./ repmat(range, size(testX, 1), 1);

Objections_1 = importdata('y_1.txt');
Objections_2 = importdata('y_2.txt');
dmodel_y1 =  polyfitn(Decision_variables, Objections_1 , order_1);
dmodel_y2 =  dacefit(Decision_variables, Objections_2, @regpoly2, @corrgauss, theta, lob, upb);

%% Optimization(NSGA-II)
nVar = size(Decision_variables,2);
VarSize = [1 nVar];
VarMin = 0;
VarMax = 1; 
nObj = numel(CostFunction(unifrnd(VarMin, VarMax, VarSize),dmodel_y1,dmodel_y2));
% NSGA-II Parameters
MaxIt = 100; 
nPop = 100;
pCrossover = 0.7;
nCrossover = 2*round(pCrossover*nPop/2);
pMutation = 0.4;
nMutation = round(pMutation*nPop);
mu = 0.02;
sigma = 0.1*(VarMax-VarMin);
% Initialization
empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];
pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    
    pop(i).Cost = CostFunction(pop(i).Position,dmodel_y1,dmodel_y2);
    
end

% Non-Dominated Sorting
[pop, F] = NonDominatedSorting(pop);

% Calculate Crowding Distance
pop = CalcCrowdingDistance(pop, F);

% Sort Population
[pop, F] = SortPopulation(pop);

% NSGA-II Main Loop

for it = 1:MaxIt
    
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2
        
        i1 = randi([1 nPop]);
        p1 = pop(i1);
        
        i2 = randi([1 nPop]);
        p2 = pop(i2);
        
        [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position);
        
        popc(k, 1).Position = (popc(k, 1).Position>VarMax).*VarMax+(popc(k, 1).Position<VarMin).*VarMin...
            +((popc(k, 1).Position<=VarMax)&(popc(k, 1).Position>=VarMin)).*popc(k, 1).Position;
        popc(k, 2).Position = (popc(k, 2).Position>VarMax).*VarMax+(popc(k, 2).Position<VarMin).*VarMin...
            +((popc(k, 2).Position<=VarMax)&(popc(k, 2).Position>=VarMin)).*popc(k, 2).Position;
        
        popc(k, 1).Cost = CostFunction(popc(k, 1).Position,dmodel_y1,dmodel_y2);
        popc(k, 2).Cost = CostFunction(popc(k, 2).Position,dmodel_y1,dmodel_y2);
        
    end
    popc = popc(:);
    
    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation
        
        i = randi([1 nPop]);
        p = pop(i);
        
        popm(k).Position = Mutate(p.Position, mu, sigma);
        popm(k).Position = (popm(k).Position>VarMax).*VarMax+(popm(k).Position<VarMin).*VarMin...
            +((popm(k).Position<=VarMax)&(popm(k).Position>=VarMin)).*popm(k).Position;
        popm(k).Cost = CostFunction(popm(k).Position,dmodel_y1,dmodel_y2);
        
    end
    
    % Merge
    pop = [pop
         popc
         popm]; %#ok
     
    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    pop = SortPopulation(pop);
    
    % Truncate
    pop = pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop = CalcCrowdingDistance(pop, F);

    % Sort Population
    [pop, F] = SortPopulation(pop);
    
    % Store F1
    F1 = pop(F{1});
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(3);
    PlotCosts(F1);
    F2 = [F1.Cost]';
    pause(0.01);
    
end

%注意，求最大值，要加负号
for i=1:numel(F1)
    ParetoSet(i,:)=F1(i).Position;
    ParetoFront(i,:)=F1(i).Cost';
end

%% Knee points
PF_norm = mapminmax(ParetoFront',0,1);
PF_norm = PF_norm';
D = sum(PF_norm.^2,2);
[~,id] = min(D);
best_solution = ParetoSet(id,:);
best_solution = variables_lob+best_solution.*range;
best_PF = ParetoFront(id,:);

%% Costfunction
function Y=CostFunction(X,dmodel_y1,dmodel_y2)
Y(1)=polyvaln(dmodel_y1,X);
Y(2)=predictor(X,dmodel_y2);
Y=Y';
end

