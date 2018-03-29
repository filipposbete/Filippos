%%%%%%%%My Code%%%%%%%%%
tic;
clear all
clc
%% Loading data
load('Returns.mat')
monthlyReturns = AverageValueWeightedReturnsMonthly;
clear 'AverageValueWeightedReturnsMonthly';
monthlyReturns(monthlyReturns == -99.99 ) =nan;
monthlyReturns=monthlyReturns/100;
%Size of Industry
load('SIZE.mat');
sizeAssets = Size;
clear Size;
%Correct data for NaN elements with value -99.99
sizeAssets ( sizeAssets==-99.99 ) =NaN;
%Normalization
 for i=1:size(sizeAssets,1)
    sizeAssets(i,2:end) = sizeAssets(i,2:end)./max(sizeAssets(i,2:end));
 end
 
%Defining sizes
%Define the number of Assets and number of total months
numAssets = size(monthlyReturns,2);
numMonths = size(monthlyReturns,1 );

%% Construction the skewnessMatrix

%We construct the Skewness Matrix based on the skewness criterio 5 year
% Considering the investment day is the 60th month hence the realized
% return refers to 61th month.
n=60;
skewnessMatrix = nan(numMonths,numAssets);
for numAsset = 1 : numAssets
   for i = 60:numMonths
    skewnessMatrix(i, numAsset) = ((sum((monthlyReturns(i-1, numAsset) - nanmean(monthlyReturns, 1)).^3))/n)/(((sum((monthlyReturns(i-1, numAsset) - nanmean(monthlyReturns, 1)).^2))/n).^1.5);
   end
end
%Resize the the matrix about sizes.
sizeAssets=imresize(sizeAssets,size(skewnessMatrix)) ;

%% Summary Statistics
treatNan = isnan(skewnessMatrix);
treatNan2 = isnan(skewnessMatrix)| isnan(sizeAssets( 1:end ,:))  ;
averageValue = nanmean( nanmean(skewnessMatrix,2) );
minValue = min(min(skewnessMatrix));
%We set initial values for the vectors we want to derive
crossValue = zeros(1,numMonths);
skew25 = crossValue ;
skewMedian = crossValue ;
skew75 = crossValue ;
skew95 = crossValue ;
skewMax = crossValue ;
numAssets = crossValue ;
skewStd =  crossValue ;
skewSkew =  crossValue ;
skewKur = crossValue ;
%Estimation of the summary statistics
for i = 60: size(skewnessMatrix,1) 
      [sortedSkewnessofReturns,indexSkewness] = sort( skewnessMatrix(i, ~treatNan(i,:)),2);
      numOfAssetsOfMonth = size(sortedSkewnessofReturns,2);
      crossValue(1,i) = prctile(sortedSkewnessofReturns,5);
      skew25(1,i)= prctile(sortedSkewnessofReturns,25);
      skewMedian(1,i) = median( sortedSkewnessofReturns(1: round(0.25*numOfAssetsOfMonth) ));
      skew75(1,i) = prctile(sortedSkewnessofReturns,75);
      skew95(1,i) =  prctile(sortedSkewnessofReturns,95);
      skewMax(1,i) = max( sortedSkewnessofReturns ); 
      numAssets(1,i) = numOfAssetsOfMonth;
      skewStd(1,i)= std( sortedSkewnessofReturns );
      skewSkew(1,i) = skewness(sortedSkewnessofReturns);
      skewKur(1,i) = kurtosis(sortedSkewnessofReturns);
      
end

min5 = mean(crossValue)
min25 = mean(skew25)
medianV = mean(skewMedian) 
min75 = mean(skew75)
min95 = mean( skew95)
maxV = mean(skewMax)
numberOfA = sum(numAssets)
stDeviationA = mean(skewStd)
skewnessA = mean(skewSkew)
kurtosisA = mean(skewKur)

%% Correlation between Momentum and Skewness
load('momentumMatrixCriterion.mat')
%Elimination of the column dates from the matrix of momentum
momentumMatrixCriterion = momentumMatrixCriterion(:,(2:end));
%We transform the matrix of skewness in order to have the same rows with
%the momentum matrix.
crossMeanMom = nanmean(momentumMatrixCriterion,2);
crossMean = nanmean(skewnessMatrix,2);
pearsonCorrelation= corrcoef(crossMeanMom,crossMean,'rows','pairwise');
spearmanCorr = spear(crossMeanMom,crossMean);

%% Creation of a sorted skewness portfolio
%Code for treatin nans
treatNanInSortCriterion = isnan(skewnessMatrix);
treatNanInSortCriterionOrSize = isnan(skewnessMatrix)| isnan(sizeAssets);
[sortedMatrixofReturns,indexMatrix] = sort( skewnessMatrix,2);
newskew =zeros(10,size(skewnessMatrix,1));
for i=60:size(skewnessMatrix,1)
    %[sorted,index] = sort(skewnessMatrix(i,~treatNan(i,:)),2);
    [sortedMatrixofReturns,indexMatrix] = sort( skewnessMatrix(i, ~treatNanInSortCriterion(i,:)),2);
    %numOfAssetsOfMonth = size(sorted,2);
    numOfAssetsOfMonth = size(sortedMatrixofReturns,2);
    skwMat(1,i) = mean(sortedMatrixofReturns(1:round(0.1*numOfAssetsOfMonth)));
    skwMat(2,i) = mean(sortedMatrixofReturns(1:round(0.1*numOfAssetsOfMonth)+1 : round(0.2*numOfAssetsOfMonth)));
    skwMat(3,i) = mean(sortedMatrixofReturns(1:round(0.2*numOfAssetsOfMonth)+1 : round(0.3*numOfAssetsOfMonth)));
    skwMat(4,i) = mean(sortedMatrixofReturns(1:round(0.3*numOfAssetsOfMonth)+1 : round(0.4*numOfAssetsOfMonth)));
    skwMat(5,i) = mean(sortedMatrixofReturns(1:round(0.4*numOfAssetsOfMonth)+1 : round(0.5*numOfAssetsOfMonth)));
    skwMat(6,i) = mean(sortedMatrixofReturns(1:round(0.5*numOfAssetsOfMonth)+1 : round(0.6*numOfAssetsOfMonth)));
    skwMat(7,i) = mean(sortedMatrixofReturns(1:round(0.6*numOfAssetsOfMonth)+1 : round(0.7*numOfAssetsOfMonth)));
    skwMat(8,i) = mean(sortedMatrixofReturns(1:round(0.7*numOfAssetsOfMonth)+1 : round(0.8*numOfAssetsOfMonth)));
    skwMat(9,i) = mean(sortedMatrixofReturns(1:round(0.8*numOfAssetsOfMonth)+1 : round(0.9*numOfAssetsOfMonth)));
    skwMat(10,i) = mean(sortedMatrixofReturns(1:round(0.9*numOfAssetsOfMonth)+1 :end));
end 
theSkewnessAverage=mean(skwMat,2);
%% %% Creation of an Equally and a Value Weighted Portfolio Returns
equallywReturns=nan(numMonths,1);
valuewReturns=nan(numMonths,1);
for iInvestmentDate=60:size(skewnessMatrix,1)-1
    %iInvestmentDate=60;
    sortCriterion=skewnessMatrix(iInvestmentDate,:);
    sortCriterion(isnan(sizeAssets(iInvestmentDate+1,:)))=nan;
    
    [sortedMatrixofReturns,indexSorted] = sort(sortCriterion ,2);
      
    numOfAssetsOfMonth=sum(~isnan(sortedMatrixofReturns),2);
    numInDecileThisIteration=round(numOfAssetsOfMonth/10);
    
    indexLongThisIteration=indexSorted(...
          numOfAssetsOfMonth-numInDecileThisIteration+1:numOfAssetsOfMonth);
      
    indexShortThisIteration=indexSorted(1:numInDecileThisIteration);
     
    equallywReturns(iInvestmentDate)...
        =nanmean(monthlyReturns(iInvestmentDate+1,indexLongThisIteration))-...
         nanmean(monthlyReturns(iInvestmentDate+1,indexShortThisIteration));      

    

     valueWeightLong= sizeAssets(iInvestmentDate,indexLongThisIteration)/sum(sizeAssets(iInvestmentDate,indexLongThisIteration));
     valueWeightShort= sizeAssets(iInvestmentDate,indexShortThisIteration)/sum(sizeAssets(iInvestmentDate,indexShortThisIteration));
     
    valuewReturns(iInvestmentDate)...
        =sum(monthlyReturns(iInvestmentDate+1,indexLongThisIteration).* valueWeightLong)-...
         sum(monthlyReturns(iInvestmentDate+1,indexShortThisIteration).*valueWeightShort);      
                                   
end

ewRwithNoNan=equallywReturns;ewRwithNoNan(isnan(ewRwithNoNan))=0;
vwRwithNoNan=valuewReturns;vwRwithNoNan(isnan(vwRwithNoNan))=0;

%% Regression for CAPM and FF models
load ('MKT.mat');
load ('SMB.mat');
load ('HML.mat');
load ('RF.mat');
%%%Correction of the size for the matrices above in order to have the same
%%%size with returns of an equally and a value weighted portfolio
minLength = min(length(equallywReturns), length(MKT));
% Removing any extra elements from the longer matrix
equallywReturns= equallywReturns(1:minLength);
MKT= MKT(1:minLength);
HML = HML(1:minLength);
SMB = SMB(1:minLength);
RF = RF(1:minLength);
%The CAPM Regression for an Equally Weighted Portfolio
[covarEqCapm,NwstdEqCapm, CapmEqBeta ] =  hac(MKT, ewRwithNoNan );
tstatEqCapm = CapmEqBeta./NwstdEqCapm
%The FF Regression for an Equally weighted Portfolio
regressors = [ MKT, SMB, HML  ];
[covarEqFAMA,NwstdEqFF, CoeffEqFF ] =  hac( regressors, ewRwithNoNan );
famaTstatisticsE = CoeffEqFF./NwstdEqFF;


%%Figures 
dm=load('DATESIZE')

figure (1)
plot(datenum(num2str(dm.Date2(1:1098)),'yyyymm'),cumprod(1+ewRwithNoNan));hold on
plot(datenum(num2str(dm.Date2(1:1098)),'yyyymm'),cumprod(1+vwRwithNoNan))
datetick

figure(2)
cumReturns = cumprod((ewRwithNoNan(60,:)/10 )+1)-1 ;
 
 
 excessReturns = ewRwithNoNan(60,:) - transpose(RF);
 cumReturns = cumprod( (excessReturns/100) +1  )-1 ;
 cumLogReturns = cumsum( log( 1+ excessReturns) ) ;
 
 cumTime =1: size(cumReturns,2);
 
 %Setting Time
 startingDate = datetime('30-June-1927');
 endingDate = datetime(2017,12,31);
 matrixOfDates = linspace( startingDate, endingDate, size(cumReturns,2)  );
 
 hold on
 
 yyaxis left
 plot(matrixOfDates,cumReturns);  
 ylabel('Cumulative Skewness%')
 
 
 % plot(cumTime,cumLogReturns);  
 yyaxis right
 plot(matrixOfDates,real(cumLogReturns));  
 legend( {'Compounded excess return' , 'Cumulative log excess return'   } );
 
 xlabel('Date')
 
hold off