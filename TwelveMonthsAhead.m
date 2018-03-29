%% %%Returns Ahead
treatNanInSortCriterion = isnan(skewnessMatrix);
treatNanInSortCriterionOrSize = isnan(skewnessMatrix)| isnan(sizeAssets)  ;

%Define Number of Categories for each Criterion
numOfCategoriesCriterion1 = 10;

%Initialization of Matrices
averageEwR = zeros(12,1);
averageVwR = averageEwR;
ValuesForTable = zeros(12,10);

for imonth = 1:12
    %Initialization of Return Matrices
    ewR=nan(numMonths,2);
    vwR=nan(numMonths,2);    
    
for i=60:numMonths - imonth
 
    sortCriterion=skewnessMatrix(i,2:end);
    sortCriterion(isnan(monthlyReturns(i+1,:)))=nan;
    [sortedMatrixofReturns,indexSorted] = sort(sortCriterion ,2);     
    numOfAssetsOfMonth=sum(~isnan(sortedMatrixofReturns),2);
    numInDecileCriterion1=round(numOfAssetsOfMonth / numOfCategoriesCriterion1);
    
    isLong=1 + indexSorted(...
          numOfAssetsOfMonth-numInDecileCriterion1+1:numOfAssetsOfMonth);      
    isShort=1 + indexSorted(1:numInDecileCriterion1);
     ewR(i+1,2)...
        =nanmean(monthlyReturns(i+imonth,isLong))-...
         nanmean(monthlyReturns(i+imonth,isShort));           
     weightLong= sizeAssets(i,isLong)/sum(sizeAssets(i,isLong));
     weightShort= sizeAssets(i,isShort)/sum(sizeAssets(i,isShort));     
    vwR(i+1,2)...
        =sum(monthlyReturns(i+imonth,isLong).* weightLong)-...
         sum(monthlyReturns(i+imonth,isShort).*weightShort);      
end
    averageEwR(imonth) = nanmean(ewR(:,2))*100;
    averageVwR(imonth) = nanmean(vwR(:,2))*100;    

    %Regression Part
    ewRExcessReturns = ewR(:,2) - (RF/100);
    vwRExcessReturns = vwR(:,2) - (RF/100);
 
    treatNaN = isnan(ewRExcessReturns);


    %Equally Weighted
    %CAPM
    [covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( MKT,ewRExcessReturns );
    capmTstatisticsE = capmCoeffE(1)/capmNwErrorsE(1);
 
    %FAMA - MACBETH
    famaIndependentVar = [MKT, SMB, HML];
    [covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar,ewRExcessReturns );
    famaTstatisticsE = famaCoeffE(1)/famaNwErrorsE(1);
 
    %Value Weighted
    %CAPM
    treatNaN = isnan(vwRExcessReturns);
    [covarCapm,capmNwErrorsV, capmCoeffV ] =  hac(MKT, vwRExcessReturns );
    capmTstatisticsV = capmCoeffV(1)/capmNwErrorsV(1);
 
    %FAMA - MACBETH
    famaIndependentVar = [ MKT, SMB, HML  ];
    [covarFAMA,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar,vwRExcessReturns );
    famaTstatisticsV = famaCoeffV(1)/famaNwErrorsV(1);
 
       
    ValuesForTable(imonth,1) = capmCoeffE(1);
    ValuesForTable(imonth,2) = capmTstatisticsE;
    ValuesForTable(imonth,3) = famaCoeffE(1);
    ValuesForTable(imonth,4) = famaTstatisticsE;
    ValuesForTable(imonth,5) = capmCoeffV(1);
    ValuesForTable(imonth,6) = capmTstatisticsV;
    ValuesForTable(imonth,7) = famaCoeffV(1);
    ValuesForTable(imonth,8) = famaTstatisticsV;
    ValuesForTable(imonth,9) = averageEwR(imonth);
    ValuesForTable(imonth,10) = averageVwR(imonth);
end