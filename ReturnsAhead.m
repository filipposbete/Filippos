averageEwR = zeros(12,1);
averageVwR = averageEwR;
statisticsTable = zeros(12,8);

for imonthsAhead = 1:12
    %Initialization of Return Matrices
    for iInvestDate = 60:numMonths-imonthsAhead
    [sortSkewness,ixSkew] = sort(skewnessMatrix,2);
    numAssetsMonth = sum(~isnan(skewnessMatrix(iInvestDate,:)));
    numDeciles = round(numAssetsMonth/numCategories);
    
    position1Decile = ixSkew(iInvestDate,1:numDeciles);
    position10decile =  ixSkew(iInvestDate,9*numDeciles+1 :numAssetsMonth);
    
    equallyMonthsAhead(iInvestDate+1,imonthsAhead)=...
         - nanmean(monthlyReturns(iInvestDate+imonthsAhead,position1Decile))+...
          nanmean(monthlyReturns(iInvestDate+imonthsAhead,position10decile)); 
      
      weight1decile = sizeAssets(iInvestDate,position1Decile)/sum(sizeAssets(iInvestDate,position1Decile));
    weight10decile = sizeAssets(iInvestDate,position10decile)/sum(sizeAssets(iInvestDate,position10decile));
    
    valueMonthsAhead(iInvestDate+1,imonthsAhead)...
        = -nansum(monthlyReturns(iInvestDate+imonthsAhead,position1Decile).* weight1decile)+...
          nansum(monthlyReturns(iInvestDate+imonthsAhead,position10decile).*weight10decile); 
    end
 theMeanEquallyMonthsAhead =  nanmean(equallyMonthsAhead,1);
theMeanValueMonthsAhead =  mean(valueMonthsAhead,1);

%FF Regresion Analysis for equally weighted
 ffIndependentVar = [ MKT, SMB, HML];
 [~,famaNwErrorsEqAhead, famaCoeffEqAhead ] =  hac( ffIndependentVar,equallyMonthsAhead(:,imonthsAhead ));
 famaTstatisticsEqAhead = famaCoeffEqAhead(1)/famaNwErrorsEqAhead(1);
 
 onesRegression = ones(numMonths,1);
 [~,onesNewWestAheadEqually, onesCoefficientEquallyAhead] = hac(onesRegression,...
                    equallyMonthsAhead(:,imonthsAhead)*100,'intercept',false);
 onesTstatisticsEqAhead = onesCoefficientEquallyAhead(1)/onesNewWestAheadEqually(1);  
 
 % FF regression Analysis for value weighted
 [~,onesNewWestAheadValue, onesCoefficientValueAhead] =  hac(onesRegression,...
                    valueMonthsAhead(:,imonthsAhead)*100,'intercept',false);
 onesTstatisticsVlAhead = onesCoefficientValueAhead(1)/onesNewWestAheadValue(1); 
 
 [~,famaNwErrorsAheadVl, famaCoeffvalueAhead ] =  hac( ffIndependentVar,valueMonthsAhead(:,imonthsAhead ));
 famaTstatisticsVlAhead = famaCoeffvalueAhead(1)/famaNwErrorsAheadVl(1);
 
 
 statisticsTable(imonthsAhead,1) = theMeanEquallyMonthsAhead(imonthsAhead)*100;
 statisticsTable(imonthsAhead,2) = onesTstatisticsEqAhead;
 statisticsTable(imonthsAhead,3) = famaCoeffEqAhead(1)*100;
 statisticsTable(imonthsAhead,4) = famaTstatisticsEqAhead;
 statisticsTable(imonthsAhead,5) = theMeanValueMonthsAhead(imonthsAhead)*100;
 statisticsTable(imonthsAhead,6) = onesTstatisticsVlAhead;
 statisticsTable(imonthsAhead,7) = famaCoeffvalueAhead(1)*100;
 statisticsTable(imonthsAhead,8) = famaTstatisticsVlAhead;
end

xlswrite('Returns Ahead Results Equally', statisticsTable(:,1:4));   
xlswrite('Returns Ahead Results Value', statisticsTable(:,5:8));  
