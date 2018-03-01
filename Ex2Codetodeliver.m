 %Loading my data
load('DATAFinanceExercise2.mat')
load('DATAIndustriesNames.mat');

%We rename the Matrix AverageValueWeightedReturnsMontly to 'MonthlyReturns'
MonthlyReturns = AverageValueWeightedReturnsMonthly;

%We name the Mean Matrix = Mean and VarCovar = Variance-Covariance Matrix
Mean = nanmean(MonthlyReturns);
VarCovar = cov(MonthlyReturns, 'partialrows');

%Efficient Frontier and the optimal Sharpe ratio
Portf = Portfolio('AssetMean', Mean, 'AssetCovar', VarCovar, 'LowerBound', 0, 'UpperBound', 1, 'Budget', 1);
plotFrontier(Portf);
figure;
plotFrontier(Portf);
Weights = estimateMaxSharpeRatio(Portf);
[risk, ret] = estimatePortMoments(Portf, Weights);
hold on
plot(risk, ret, '*k');
title('Efficient Frontier and Sharpe Ratio');
SumOfWeights = sum(Weights);

%Presentation of weights, using pie.
TransposeWeights = reshape(Weights, 1, 49);
figure;
pie(TransposeWeights);
legend('Agric','Food','Soda','Beer','Smoke','Toys','Fun','Books','Hshld',...
       'Clths','Hlth','MedEq','Drugs','Chems','Rubbr','Txtls','BldMt','Cnstr','Steel',...
       'FabPr','Mach','ElcEq','Autos','Aero','Ships','Guns','Gold','Mines','Coal','Oil','Util',...
       'Telcm','PerSv','BusSv','Hardw','Softw','Chips','LabEq','Paper','Boxes',...
       'Trans','Whlsl','Rtail','Meals','Banks','Insur','RIEst','Fin','Other');



