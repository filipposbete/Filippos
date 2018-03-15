%Loading Data about the Size, Returns and book to market ratio
load('Returns.mat')
Returns = AverageValueWeightedReturnsMonthly;
clear AverageValueWeightedReturnsMonthly;
Returns(Returns==-99.99) = NaN;
load('SIZE.mat');
sizeIndustry = Size;
clear Size;
sizeIndustry(sizeIndustry == -99.99) = NaN;
load('booktomarket.mat');

%We set the No of assets and the No periods
noMonths = size(Returns, 1);
noAssets = size(Returns, 2);

%We construct the Skewness Matrix 
n=60;
skewnessMatrix = zeros(size(Returns, 1),size(Returns, 2));
for numAsset = 1:size(Returns, 2)
   for i = 1:size(Returns,1)
    skewnessMatrix(i, numAsset) = ((sum((Returns(i, numAsset) - nanmean(Returns, 1)).^3))/n)/(((sum((Returns(i, numAsset) - nanmean(Returns, 1)).^2))/n).^1.5);
   end
end

%Calculation of Statistics - Panel A
crossMean = nanmean(skewnessMatrix,2);
skwMean = nanmean(crossMean,1);
skwStd = std(crossMean,1);
skwKyrtosis = kurtosis(crossMean,1);
skwkewness = skewness(crossMean,1);
skwMedian = median(crossMean,1);
skwMax = max(crossMean);
skwMin = min(crossMean);
skw5percent = prctile(crossMean,0.05);
skw25percent = prctile(crossMean,0.25);
skw75percent = prctile(crossMean,0.75);
skw95percent = prctile(crossMean,0.95);

%Correlation skewness with momentum - Panel B
load('momentumMatrix109849.mat')
%We transform the matrix of skewness in order to have the same rows with
%the momentum matrix.
momentumMatrix= imresize(momentumMatrix,size(skewnessMatrix)) ;
crossMeanMom = nanmean(momentumMatrix,2);

spearmanCorr = corr(crossMeanMom,crossMean,'Type','Spearman');
pearsonCorr = corr(crossMeanMom,crossMean,'Type','Pearson');

%A shorted Skewness portfolio
treatNan = isnan(skewnessMatrix)|isnan(Returns);
skwMat = zeros(10,size(skewnessMatrix,1));
for i=1:size(skewnessMatrix,1)
    [sorted,index] = sort(skewnessMatrix(i,~treatNan(i,:)),2);
    nofAssetMo = size(sorted,2);
    skwMat(1,i) = mean(sorted(1:round(0.1*nofAssetMo)));
    skwMat(2,i) = mean(sorted(1:round(0.1*nofAssetMo)+1 : round(0.2*nofAssetMo)));
    skwMat(3,i) = mean(sorted(1:round(0.2*nofAssetMo)+1 : round(0.3*nofAssetMo)));
    skwMat(4,i) = mean(sorted(1:round(0.3*nofAssetMo)+1 : round(0.4*nofAssetMo)));
    skwMat(5,i) = mean(sorted(1:round(0.4*nofAssetMo)+1 : round(0.5*nofAssetMo)));
    skwMat(6,i) = mean(sorted(1:round(0.5*nofAssetMo)+1 : round(0.6*nofAssetMo)));
    skwMat(7,i) = mean(sorted(1:round(0.6*nofAssetMo)+1 : round(0.7*nofAssetMo)));
    skwMat(8,i) = mean(sorted(1:round(0.7*nofAssetMo)+1 : round(0.8*nofAssetMo)));
    skwMat(9,i) = mean(sorted(1:round(0.8*nofAssetMo)+1 : round(0.9*nofAssetMo)));
    skwMat(10,i) = mean(sorted(1:round(0.9*nofAssetMo)+1 :end));
end 
theAverage=mean(skwMat,2);

%Construction of Portfolio based on skewness strategy.
%Sort of Skewness Matrix i.e. sorting the returns of each t according to
%the skewness 
[sorted,index] = sort(skewnessMatrix,2);
for i=1:size(skewnessMatrix,1)  
   for j=1:size(skewnessMatrix,2)
       Returns(i,j) = Returns(i,index(i,j));
   end
end
treatNan = isnan(Returns);
%Estimation of the mean return in each decile every month
equalWeighted = zeros(1098,10);
for i=1:size(Returns,1) 
    numberOfAsset =  size(Returns,2) - sum(treatNan(i,:));
    NumofAssetInDecileRound = round( 0.1 * numberOfAsset);
    equalWeighted(i,1) = mean(Returns(i,1:round(0.1*numberOfAsset)));    
    t=0.1;
    l=0.2;
    for j=2:9
        equalWeighted(i,j) = mean(Returns(i,round(t*numberOfAsset)+1 : round(l*numberOfAsset)));
        t=t+0.1; 
        l=l+0.1;    
    end
    equalWeighted(i,10) = mean(Returns(i,round(0.9*numberOfAsset)+1 : numberOfAsset));
end
%Hence, we have the returns for each decile
for i=1:10
portMeanDecile(1,i) = mean(equalWeighted(:,i));
end
% Montlhy Returns of equally balanced portfolio using skewness
SkewnessMonthlyReturn = equalWeighted(:,1) - equalWeighted(:,10);


%Loading data for the Fama - French Model
load('MKT');
load('HML');
load('SMB');
sizeIndustry=imresize(sizeIndustry,size(skewnessMatrix)) ;
%bookToMarketMonthly=imresize(bookToMarketMonthly,size(skewnessMatrix)) ;
MKT= imresize(MKT,size(SkewnessMonthlyReturn)) ;
HML= imresize(HML,size(SkewnessMonthlyReturn)) ;
SMB= imresize(SMB,size(SkewnessMonthlyReturn)) ;

% %CAPM Model 
regStatsCAPM = regstats(SkewnessMonthlyReturn,MKT, 'linear');
BetaCAPM = regStatsCAPM.tstat.beta;
ResidaulCAPM = regStatsCAPM.r;
newSTDCAPM = nwse(ResidaulCAPM,MKT);
tstCAPM = BetaCAPM ./ newSTDCAPM;

%Fama-French Model
regressors = [MKT, SMB, HML];
regStatsFFF = regstats(SkewnessMonthlyReturns,regressors, 'linear');
BetaFFF = regStatsFFF .tstat.beta;
ResidaulFFF = regStatsFFF .r;
newstandarerrorsEq = nwse(ResidaulFFF,Y);
tstatFFFEq = BetaFFF ./ newstandarerrorsEq;


