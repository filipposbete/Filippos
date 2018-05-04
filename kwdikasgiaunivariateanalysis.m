%We sort the skewness Matrix
[sortSkewness,ixSkew] = sort(skewnessMatrix,2);

%We construct some matrices
realizedReturnsEqually = nan(numMonths,11);
realizedReturnsValue = nan(numMonths,11);
numCategories = 10;
for iInvestDate = 60: numMonths-1
    [sortSkewness,ixSkew] = sort(skewnessMatrix,2);
    numAssetsMonth = sum(~isnan(skewnessMatrix(iInvestDate,:)));
    numDeciles = round(numAssetsMonth/numCategories);
%     %For the 1st decile
%     longPosition1decile = ixSkew(iInvestDate,1 : numDeciles);
%     realizedReturnsEqually(iInvestDate+1,1) = nanmean((monthlyReturns(iInvestDate+1,longPosition1decile)));
%    
    %For the 1st to 9th decile
    %The 1st decile refers to decile that we go adopt a long position
    %Based on the negative skewness Strategy
    temp = 0;
     for iDecile=0:1:8
         temp = temp +1;
        position = ixSkew(iInvestDate,(iDecile*numDeciles)+1 : (iDecile+1)*numDeciles);
        realizedReturnsEqually(iInvestDate+1,temp)= nanmean((monthlyReturns(iInvestDate+1,position)));
     end

    %For the 10th decile
    position10decile =  ixSkew(iInvestDate,9*numDeciles + 1 : numAssetsMonth);
    realizedReturnsEqually(iInvestDate+1,10) = nanmean((monthlyReturns(iInvestDate+1,position10decile)));
    
    
end


%Now we construct the code for the value weighted portfolio
realizedReturnsValue = nan(numMonths,11);
numCategories = 10;
for iInvestDate = 60 : numMonths-1
  
    [sortSkewness,ixSkew] = sort(skewnessMatrix,2);
    numAssetsMonth = sum(~isnan(skewnessMatrix(iInvestDate,:)));
    numDeciles = round(numAssetsMonth/numCategories);
%     %For the 1st decile
%     longPosition1decile = ixSkew(iInvestDate,1:numDeciles);
%     weightDecile1decile = sizeAssets(iInvestDate,longPosition1decile)/...
%             nansum(sizeAssets(iInvestDate,longPosition1decile)); 
%     realizedReturnsValue(iInvestDate+1,1)= nansum(monthlyReturns(iInvestDate+1,longPosition1decile)...
%                                                 .*weightDecile1decile);
    temp = 0;
     for iDecile=0:1:8
         temp = temp +1;
        position = ixSkew(iInvestDate,(iDecile*numDeciles)+1 : (iDecile+1)*numDeciles);
        
        weightDecile = sizeAssets(iInvestDate,position)/...
            nansum(sizeAssets(iInvestDate,position)); 
       
        realizedReturnsValue(iInvestDate+1,temp)= nansum(monthlyReturns(iInvestDate+1,position)...
                                                .*weightDecile);
    end
     
     %For the 10th decile
    position10decile = ixSkew(iInvestDate,9*numDeciles +1 :numAssetsMonth);
    
    weightShort10decile= sizeAssets(iInvestDate,position10decile)/...
                    nansum(sizeAssets(iInvestDate,position10decile));
    
    realizedReturnsValue(iInvestDate+1,10)= nansum(monthlyReturns(iInvestDate+1,position10decile)...
                                            .*weightShort10decile);
   
%    %Value -Weighted Returns                                     
%    realizedReturnsValue(iInvestDate+1,11) =  realizedReturnsValue(iInvestDate+1,1)...
%                                     - realizedReturnsValue(iInvestDate+1,10);
%    
end


%Now the averages in deciles ans the returns 10-1
%table 2 Panel A
for i=1:10
    equallyMeanDecile(1,i) = nanmean(realizedReturnsEqually(:,i));
end
%In percentage 
equallyMeanDecile = equallyMeanDecile*100

realizedReturnsEqually(:,11)=-realizedReturnsEqually(:,1)+realizedReturnsEqually(:,10);
%In percentage
returnsEqually = realizedReturnsEqually(:,11)*100;

realizedReturnsEquallyNoNan=realizedReturnsEqually(:,11);realizedReturnsEquallyNoNan...
            (isnan(realizedReturnsEquallyNoNan))=0;

for i=1:10
valueMeanDecile(1,i) = nanmean(realizedReturnsValue(:,i));
end

%In percentage
valueMeanDecile = valueMeanDecile*100
realizedReturnsValue(:,11) = - realizedReturnsValue(:,1) +realizedReturnsValue(:,10); 
%In percentage
returnsValue = realizedReturnsValue(:,11)*100;

realizedReturnsValueNoNan=realizedReturnsValue(:,11);realizedReturnsValueNoNan...
            (isnan(realizedReturnsValueNoNan))=0;
        
avergeRealizedReturnsEqually = (mean(realizedReturnsEquallyNoNan))*100;
averageRealizedReturnsValue = (mean(realizedReturnsValueNoNan))*100;

dm=load('DATESIZE')        
 figure
plot(datenum(num2str(dm.Date2(60:1098)),'yyyymm'),cumprod(1+realizedReturnsValueNoNan(60:1098)),':');
hold on
plot(datenum(num2str(dm.Date2(60:1098)),'yyyymm'),cumprod(1+realizedReturnsEquallyNoNan(60:1098)),'-.');
datetick


%Regresion Analysis for the returns of the equally...
%and value weigthed portfolio

%Loading 
load('MKT.mat');
MKT = MKT(1:1098,:);
load('SMB.mat');
SMB = SMB(1:1098,:);
load('HML.mat');
HML = HML(1:1098,:);

treatNaNeq = isnan(returnsEqually);
treatNaNval = isnan(returnsValue);
%CAPM Regression Analysis - Equally
[~, nwStdEqCapm, betasEqCapm ] =  hac(MKT(~treatNaNeq), returnsEqually(~treatNaNeq) );
tstEqCapmA = betasEqCapm(1)./nwStdEqCapm(1);
%CAPM Regression Analysis - Value
[~, nwStdVwCapm, betasVwCapm ] =  hac(MKT(~treatNaNval), returnsValue(~treatNaNval));
tstVwCapmA = betasVwCapm(1)./nwStdVwCapm(1);
capmAbetaEqually = betasEqCapm (1,:);
capmAbetaValue = betasVwCapm(1,:);

%FF Regression Analysis - Equally
independentVariablesEq = [ MKT(~treatNaNeq), SMB(~treatNaNeq), HML(~treatNaNeq)  ];

[~,famaNwStdEq, betasEqFf ] =  hac( independentVariablesEq,returnsEqually(~treatNaNeq));
tstFfEq = betasEqFf(1)./famaNwStdEq(1);

%FF regression Analysis - Value
independentVariablesValue = [ MKT(~treatNaNval), SMB(~treatNaNval), HML(~treatNaNval)  ];
[~,famaNwStdVw, betasVwFf ] =  hac(independentVariablesValue,returnsValue(~treatNaNval));
tstFfVl= betasVwFf(1)./famaNwStdVw(1);

FfAbetaEqually = betasEqFf(1,:);
FfAbetaValue = betasVwFf (1,:);

 onesRegression = ones(numMonths,1);
%Ones Regression for equally weighted portfolio
[~,onesNewWestErrorEqually, onesCoefficientEqually] = hac(onesRegression(~treatNaNeq),...
                    returnsEqually(~treatNaNeq)*100,'intercept',false);
 onesTstatisticsEq = onesCoefficientEqually(1)/onesNewWestErrorEqually(1);
 
 %Ones regresion for the value weighted portfolio
[~,onesNewWestErrorValue, onesCoefficientValue] =  hac(onesRegression(~treatNaNval),...
                    returnsValue(~treatNaNval)*100,'intercept',false);
 onesTstatisticsVl = onesCoefficientValue(1)/onesNewWestErrorValue(1); 
 
%
avergeRealizedReturnsEqually = (nanmean(realizedReturnsEquallyNoNan))*100;
averageRealizedReturnsValue = (nanmean(realizedReturnsValueNoNan))*100;

%Printing results for the equally weighted portfolio analysis
matrixExcelEqually = [avergeRealizedReturnsEqually ; onesTstatisticsEq ; capmAbetaEqually;...
     tstEqCapmA ; FfAbetaEqually ; tstFfEq; ];
xlswrite('Univariate Portofolio Equally',matrixExcelEqually);

%Printing results for the value weighted portfolio analysis
matrixExcelValue = [averageRealizedReturnsValue ; onesTstatisticsVl ; capmAbetaValue ;...
        tstVwCapmA ; FfAbetaValue ; tstFfVl ;];
xlswrite('Univariate Portofolio Value',matrixExcelValue);   
   
xlswrite('Returns for Equally Portfolio', returnsEqually)