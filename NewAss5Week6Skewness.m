tic;
%Loading Data about the Size, Returns and book to market ratio
clear all;
load('Returns.mat')
monthlyReturns = AverageValueWeightedReturnsMonthly;
clear AverageValueWeightedReturnsMonthly;
monthlyReturns(monthlyReturns==-99.99) = NaN;
load('SIZE.mat');
sizeIndustry = Size;
clear Size;
sizeIndustry(sizeIndustry == -99.99) = NaN;
load('booktomarket.mat');
%Creating monthly book to market matrix
bookToMarketMonthly = zeros( 12*size(booktomarket,1), size(monthlyReturns,2 ) );
temp=1;
positionInBookToMarketMatrix =1;
for i=1:12*size(booktomarket,1)
    if(temp==13)
        temp=1;
        positionInBookToMarketMatrix =positionInBookToMarketMatrix +1;
    end
    bookToMarketMonthly(i,:)=booktomarket(positionInBookToMarketMatrix,:) ;
    temp=temp+1;
end
%We set the No of assets and the No periods
noMonths = size(monthlyReturns, 1);
noAssets = size(monthlyReturns, 2);

%We construct the Skewness Matrix 
n=60;
skewnessMatrix = zeros(size(monthlyReturns, 1),size(monthlyReturns, 2));
for numAsset = 1:size(monthlyReturns, 2)
   for i = 1:size(monthlyReturns,1)
    skewnessMatrix(i, numAsset) = ((sum((monthlyReturns(i, numAsset) - nanmean(monthlyReturns, 1)).^3))/n)/(((sum((monthlyReturns(i, numAsset) - nanmean(monthlyReturns, 1)).^2))/n).^1.5);
   end
end
load('size.mat');  
sizeAssets = Size(1:1098,:);
 %Correct data for NaN elements with value -99.99
 sizeAssets ( sizeAssets==-99.99 ) =NaN;
 for i=1:size(sizeAssets,1)
    sizeAssets(i,:) = sizeAssets(i,:)./max(sizeAssets(i,:));
 end
treatNan = isnan(skewnessMatrix);
treatNan2 = isnan(skewnessMatrix)| isnan(sizeAssets( 1:end ,:))  ;
averageValue = mean( nanmean(skewnessMatrix,2) );
minValue = min(min(skewnessMatrix));
%We set initial values for the vectors we want to derive
crossValue = zeros(1,noMonths);
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
for i = 1: size(skewnessMatrix,1) 
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

%Correlation skewness with momentum - Panel B
load('momentumMatrix109849.mat')
%We transform the matrix of skewness in order to have the same rows with
%the momentum matrix.
momentumMatrix= imresize(momentumMatrix,size(skewnessMatrix)) ;
crossMeanMom = nanmean(momentumMatrix,2);
crossMean = nanmean(skewnessMatrix,2)
spearmanCorr = corr(crossMeanMom,crossMean,'Type','Spearman');
pearsonCorr = corr(crossMeanMom,crossMean,'Type','Pearson');




%Sorted Portfolio Analysis
%A shorted Skewness portfolio
%treatNan = isnan(skewnessMatrix)|isnan(monthlyReturns);
%skwMat = zeros(10,size(skewnessMatrix,1));
treatNan = isnan(skewnessMatrix);
[sortedSkewnessReturns,indexSkewness] = sort(skewnessMatrix,2);
newSkew =zeros(10,size(skewnessMatrix,1));
for i=1:size(skewnessMatrix,1)
    [sortedSkewnessReturns,indexSkewness] = sort( skewnessMatrix(i, ~treatNan(i,:)),2);   
    numOfAssetsOfMonth = size(sortedSkewnessReturns,2);
    skwMat(1,i) = mean(sortedSkewnessReturns(1:round(0.1*numOfAssetsOfMonth)));
    skwMat(2,i) = mean(sortedSkewnessReturns(1:round(0.1*numOfAssetsOfMonth)+1 : round(0.2*numOfAssetsOfMonth)));
    skwMat(3,i) = mean(sortedSkewnessReturns(1:round(0.2*numOfAssetsOfMonth)+1 : round(0.3*numOfAssetsOfMonth)));
    skwMat(4,i) = mean(sortedSkewnessReturns(1:round(0.3*numOfAssetsOfMonth)+1 : round(0.4*numOfAssetsOfMonth)));
    skwMat(5,i) = mean(sortedSkewnessReturns(1:round(0.4*numOfAssetsOfMonth)+1 : round(0.5*numOfAssetsOfMonth)));
    skwMat(6,i) = mean(sortedSkewnessReturns(1:round(0.5*numOfAssetsOfMonth)+1 : round(0.6*numOfAssetsOfMonth)));
    skwMat(7,i) = mean(sortedSkewnessReturns(1:round(0.6*numOfAssetsOfMonth)+1 : round(0.7*numOfAssetsOfMonth)));
    skwMat(8,i) = mean(sortedSkewnessReturns(1:round(0.7*numOfAssetsOfMonth)+1 : round(0.8*numOfAssetsOfMonth)));
    skwMat(9,i) = mean(sortedSkewnessReturns(1:round(0.8*numOfAssetsOfMonth)+1 : round(0.9*numOfAssetsOfMonth)));
    skwMat(10,i) = mean(sortedSkewnessReturns(1:round(0.9*numOfAssetsOfMonth)+1 :end));
end 
theAverage=mean(skwMat,2);
theAverage = reshape(theAverage,1,10);
%Sorting montlhy returns according to skewness 
[sorted,index] = sort(skewnessMatrix,2);
for i=1:size(skewnessMatrix,1)  
   for j=1:size(skewnessMatrix,2)
       Returns(i,j) = skewnessMatrix(i,index(i,j));
   end
end

%The Nan values are in the last rows in each line in matrix 'Returns'
treatNan = isnan(Returns);
%Estimation the mean return in each decile every month
equallyMatrix = zeros(1098,10);
weightMatrix = zeros(1098,10);
for i=1:size(Returns,1) 
    noAssetPerMonth = size(Returns,2) - sum(treatNan(i,:));
    noAssetDecile = round( 0.1 * noAssetPerMonth);
    equallyMatrix(i,1) = mean(Returns(i,1:round(0.1*noAssetPerMonth)));    
    t=0.1;
    l=0.2;
    for j=2:9
        equallyMatrix(i,j) = mean(Returns(i,round(t*noAssetPerMonth)+1 : round(l*noAssetPerMonth)));
        t=t+0.1; 
        l=l+0.1;    
    end
    equallyMatrix(i,10) = mean(Returns(i,round(0.9*noAssetPerMonth)+1 : noAssetPerMonth));
end

%table 2 Panel A
for i=1:10
equallyMeanDecilethroughTime(1,i) = mean(equallyMatrix(:,i));
end


[sorted,index] = sort(skewnessMatrix,2);
treatNan = isnan(Returns);
weightMatrix = zeros(1098,10);
weightMatrix = Returns .* sorted ;    
for i=1:size(Returns,1) 
    noAssetPerMonth =  size(Returns,2) - sum(treatNan(i,:));
    noAssetDecile = round( 0.1 * noAssetPerMonth);
    weightMatrix(i,1) = sum(weightMatrix(i,1:round(0.1*noAssetPerMonth))) / sum(sorted(i,1:round(0.1*noAssetPerMonth)));   
    t=0.1;
    l=0.2;
    for j=2:9
        weightMatrix(i,j) = sum(weightMatrix(i,round(t*noAssetPerMonth)+1 : round(l*noAssetPerMonth))) / sum(sorted(i,round(t*noAssetPerMonth)+1 : round(l*noAssetPerMonth)));
        t=t+0.1; 
        l=l+0.1;    
    end
    weightMatrix(i,10) = sum(weightMatrix(i,round(0.9*noAssetPerMonth)+1 : noAssetPerMonth)) / sum(sorted(i,round(0.9*noAssetPerMonth)+1 : noAssetPerMonth));
end
for i=1:10
WeightMeanDecilethroughTime(1,i) = mean(weightMatrix(:,i));
WeightMeanDecilethroughTime(1,i) = WeightMeanDecilethroughTime(1,i)/100;
end
% Montlhy Returns of equally balanced and weight balanced skewness portfolio 
skewnessEquallyReturn = equallyMatrix(:,1) - equallyMatrix(:,10);
skewnessWeightReturns =  weightMatrix(:,1) -  weightMatrix(:,10);   
skewnessWeightReturns = skewnessWeightReturns/100;
%Loading data for the Fama - French Model
load('MKT');
load('HML');
load('SMB');
sizeIndustry=imresize(sizeIndustry,size(skewnessMatrix)) ;
%bookToMarketMonthly=imresize(bookToMarketMonthly,size(skewnessMatrix)) ;
MKT= imresize(MKT,size(skewnessEquallyReturn)) ;
HML= imresize(HML,size(skewnessEquallyReturn)) ;
SMB= imresize(SMB,size(skewnessEquallyReturn)) ;

% %CAPM Model - Equally Weighted
regStatsCAPM = regstats(skewnessEquallyReturn,MKT, 'linear');
BetaCAPM = regStatsCAPM.tstat.beta;
ResidualCAPM = regStatsCAPM.r;
newSTDCAPM = nwse(ResidualCAPM,MKT);
tstCAPM = BetaCAPM ./newSTDCAPM;

%Fama-French Model - Equally Weighted
regressors = [MKT, SMB, HML];
regStatsFFF = regstats(skewnessEquallyReturn,regressors, 'linear');
BetaFFF = regStatsFFF .tstat.beta;
ResidaulFFF = regStatsFFF.r;
newstandarerrorsEq = nwse(ResidaulFFF,regressors);
tstatFFFEq = BetaFFF./newstandarerrorsEq;

%CAPM Model - Value Weighted
regValueStatsCAPM = regstats(skewnessWeightReturns,MKT, 'linear');
BetaValueCAPM = regValueStatsCAPM.tstat.beta;
ResidualValueCAPM = regValueStatsCAPM.r;
newstandarerrorsValueCAPM = nwse(ResidualValueCAPM,MKT);
tstatValueCAPM = BetaValueCAPM./newstandarerrorsValueCAPM;

regressors = [MKT, SMB, HML];
%Fama French - Value Weighted
regValueStatsFFF = regstats(skewnessWeightReturns,regressors, 'linear');
BetaFFFValue = regValueStatsFFF .tstat.beta;
ResidaulFFFValue = regValueStatsFFF.r;
newstandarerrorsValueFFF = nwse(ResidaulFFFValue,MKT);
tstatValueFFF = BetaValueCAPM./newstandarerrorsValueFFF;

%Bivariate Dependent - Sort Portfolio Analysis control for size and
%momentum
    %skewnessMatrix & SizeMatrix
    [sortedSkewMatrix,ind] = sort(skewnessMatrix,2);
    for r = 1:size(skewnessMatrix,1)
       skewSIZE(r,:) = sizeIndustry(r,ind(r,:));
    end
    
    for i=1:size(skewSIZE,1)  
       for j=1:size(skewSIZE,2)
           skewSizeReturns(i,j) = monthlyReturns(i,index(i,j));
   end
end
%The Nan values are in the last rows in each line in matrix 'Returns'
treatNan = isnan(skewSizeReturns);
%Estimation the mean return in each decile every month
sequallyMatrix = zeros(1098,5);
weightMatrix = zeros(1098,5);
for i=1:size(skewSizeReturns,1) 
    noAssetPerMonth = size(skewSizeReturns,2) - sum(treatNan(i,:));
    noAssetDecile = round( 0.2 * noAssetPerMonth);
    sequallyMatrix(i,1) = mean(skewSizeReturns(i,1:round(0.2*noAssetPerMonth)));    
    t=0.1;
    l=0.2;
    for j=2:5
        sequallyMatrix(i,j) = mean(skewSizeReturns(i,round(t*noAssetPerMonth)+1 : round(l*noAssetPerMonth)));
        t=t+0.1; 
        l=l+0.1;    
    end
    sequallyMatrix(i,5) = mean(skewSizeReturns(i,round(0.8*noAssetPerMonth)+1 : noAssetPerMonth));
end

%table 2 Panel A
for i=1:5
sequallyMeanDecilethroughTime(1,i) = mean(sequallyMatrix(:,i));
end
% skewSIZE

% Value - Weighted analysis for size ans skewness
% treatNan = isnan(skewSizeReturns);
% sweightMatrix = zeros(1098,10);
% sweightMatrix = skewSizeReturns .* skewSIZE ;    
% for i=1:size(skewSizeReturns,1) 
%     noAssetPerMonth =  size(skewSizeReturns,2) - sum(treatNan(i,:));
%     noAssetDecile = round( 0.2 * noAssetPerMonth);
%     sweightMatrix(i,1) = sum(sweightMatrix(i,1:round(0.2*noAssetPerMonth))) / sum(skewSIZE(i,1:round(0.1*noAssetPerMonth)));   
%     t=0.1;
%     l=0.2;
%     for j=2:5
%         sweightMatrix(i,j) = sum(sweightMatrix(i,round(t*noAssetPerMonth)+1 : round(l*noAssetPerMonth))) / sum(skewSIZE(i,round(t*noAssetPerMonth)+1 : round(l*noAssetPerMonth)));
%         t=t+0.1; 
%         l=l+0.1;    
%     end
%     sweightMatrix(i,5) = sum(sweightMatrix(i,round(0.8*noAssetPerMonth)+1 : noAssetPerMonth)) / sum(skewSIZE(i,round(0.8*noAssetPerMonth)+1 : noAssetPerMonth));
% end
% for i=1:5
% sWeightMeanDecilethroughTime(1,i) = mean(sweightMatrix(:,i));
% sWeightMeanDecilethroughTime(1,i) = sWeightMeanDecilethroughTime(1,i)/100;
% end
% Montlhy Returns of equally balanced and weight balanced skewness portfolio 
sskewnessEquallyReturn = sequallyMatrix(:,1) - sequallyMatrix(:,5);
% sskewnessWeightReturns =  sweightMatrix(:,1) -  sweightMatrix(:,5);   
% sskewnessWeightReturns = sskewnessWeightReturns/100;

% %CAPM Model - Equally Weighted
regStatsskewSizeCAPM = regstats(sskewnessEquallyReturn,MKT, 'linear');
BetaskewSizeCAPM = regStatsskewSizeCAPM.tstat.beta;
ResidualskewSizeCAPM = regStatsskewSizeCAPM.r;
newSTDskewSizeCAPM = nwse(ResidualskewSizeCAPM,MKT);
tstskewSizeCAPM = BetaskewSizeCAPM ./newSTDskewSizeCAPM;

%Fama-French Model - Equally Weighted
regressors = [MKT, SMB, HML];
regStatsskewSizeFFF = regstats(sskewnessEquallyReturn,regressors, 'linear');
BetaskewSizeFFF = regStatsskewSizeFFF .tstat.beta;
ResidaulskewSizeFFF = regStatsskewSizeFFF.r;
newstandarerrorsskewSizeEq = nwse(ResidaulskewSizeFFF,regressors);
tstatskewSizeFFFEq = BetaskewSizeFFF./newstandarerrorsskewSizeEq;

%Sort based on skewness and book to market
%skewnessMatrix & Book to market Matrix
    [sortedBKTMatrix,ind] = sort(skewnessMatrix,2);
    for r = 1:size(bookToMarketMonthly,1)
       skewBKT(r,:) = bookToMarketMonthly(r,ind(r,:));
    end
    [sorterSkewBKT,indexBKT] = sort(skewBKT,2);
    for i=1:size(skewBKT,1)  
       for j=1:size(skewBKT,2)
           skewBKTReturns(i,j) = monthlyReturns(i,indexBKT(i,j));
   end
    end
    %The Nan values are in the last rows in each line in matrix 'Returns'
treatNan = isnan(skewBKTReturns);
%Estimation the mean return in each decile every month
BKTequallyMatrix = zeros(1098,5);
%weightMatrix = zeros(1098,5);
for i=1:size(skewBKTReturns,1) 
    noAssetPerMonth = size(skewBKTReturns,2) - sum(treatNan(i,:));
    noAssetDecile = round( 0.2 * noAssetPerMonth);
    BKTequallyMatrix(i,1) = mean(skewBKTReturns(i,1:round(0.2*noAssetPerMonth)));    
    t=0.1;
    l=0.2;
    for j=2:5
        BKTequallyMatrix(i,j) = mean(skewBKTReturns(i,round(t*noAssetPerMonth)+1 : round(l*noAssetPerMonth)));
        t=t+0.1; 
        l=l+0.1;    
    end
    BKTequallyMatrix(i,5) = mean(skewBKTReturns(i,round(0.8*noAssetPerMonth)+1 : noAssetPerMonth));
end

%table 2 Panel A
for i=1:5
BKTequallyMeanDecilethroughTime(1,i) = mean(BKTequallyMatrix(:,i));
end
BKTEquallyReturn = BKTequallyMatrix(:,1) - BKTequallyMatrix(:,5);
% sskewnessWeightReturns =  sweightMatrix(:,1) -  sweightMatrix(:,5);   
% sskewnessWeightReturns = sskewnessWeightReturns/100;
% %CAPM Model - Equally Weighted
regStatsBKTCAPM = regstats(BKTEquallyReturn,MKT, 'linear');
BetaBKTCAPM = regStatsBKTCAPM.tstat.beta;
ResidualsBKTCAPM = regStatsBKTCAPM.r;
newSTDBKTCAPM = nwse(ResidualsBKTCAPM,MKT);
tstBKTCAPM = BetaBKTCAPM ./newSTDBKTCAPM;

%Fama-French Model - Equally Weighted
regressors = [MKT, SMB, HML];
regStatsBKTFFF = regstats(BKTEquallyReturn,regressors, 'linear');
BetaBKTFFF = regStatsBKTFFF .tstat.beta;
ResidaulsBKTFFF = regStatsBKTFFF.r;
newstandarerrorsBKTEq = nwse(ResidaulsBKTFFF,regressors);
tstatBKTFFFEq = BetaBKTFFF./newstandarerrorsBKTEq;


%%Plot Cumulative Excess Returns and Log Excess Returns
load ('RF.mat')
RF = ME_BEME_RETSusingthe201801CRSPdatabase;
RF = imresize(RF,size(skewnessWeightReturns));
cumReturns = cumprod((skewnessWeightReturns(11,:)/10 )+1)-1 ;
 figure(1)
 
 excessReturns = skewnessEquallyReturn(11,:) - transpose(RF);
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
 ylabel('Cumulative MOM%')
 
 
 % plot(cumTime,cumLogReturns);  
 yyaxis right
 plot(matrixOfDates,real(cumLogReturns));  
 legend( {'Compounded excess return' , 'Cumulative log excess return'   } );
 
 xlabel('Date')
 
hold off

