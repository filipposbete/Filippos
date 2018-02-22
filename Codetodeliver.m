%My data
load('DATAFinanceExercise2.mat')
load('DATAIndustriesNames.mat');

%We rename the Matrix AverageValueWeightedReturnsMontly to A
A = AverageValueWeightedReturnsMonthly;

%We name the Mean Matrix = M and C = Variance-Covariance Matrix
M = nanmean(A);
C = cov( A ,'partialrows');

%Efficient Frontier and the optimal Sharpe ratio
p = Portfolio('AssetMean',M , 'AssetCovar', C, 'LowerBound', 0,'UpperBound',1, 'Budget',1);
plotFrontier(p);
figure;
plotFrontier(p);
portweights = estimateMaxSharpeRatio(p);
[risk, ret] = estimatePortMoments(p, portweights);
hold on
plot(risk,ret,'*k');
title('Efficient Frontier and Sharpe Ratio');
sum(portweights);

%Create a pie to present the weights.
k=reshape(portweights,1,49);
figure;
pie(k);
legend('Agric','Food','Soda','Beer','Smoke','Toys','Fun','Books','Hshld','Clths','Hlth','MedEq','Drugs','Chems','Rubbr','Txtls','BldMt','Cnstr','Steel',...
'FabPr','Mach','ElcEq','Autos','Aero','Ships','Guns','Gold','Mines','Coal','Oil','Util','Telcm','PerSv','BusSv','Hardw','Softw','Chips','LabEq','Paper','Boxes',...
'Trans','Whlsl','Rtail','Meals','Banks','Insur','RIEst','Fin','Other');



