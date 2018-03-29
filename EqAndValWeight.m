%% Creaton of an equally and value weighted portfolio in deciles
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


%treatNan = isnan(skewnessMatrix);

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
WeightMeanDecilethroughTime(1,i) = nanmean(weightMatrix(:,i));
WeightMeanDecilethroughTime(1,i) = WeightMeanDecilethroughTime(1,i)/100;
end
% Montlhy Returns of equally balanced and weight balanced skewness portfolio 
skewnessEquallyReturn = equallyMatrix(:,1) - equallyMatrix(:,10);
skewnessWeightReturns =  weightMatrix(:,1) -  weightMatrix(:,10);   
skewnessWeightReturns = skewnessWeightReturns/100;

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