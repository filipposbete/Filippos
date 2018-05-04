
%We define the sorting criteria
%Considering the size
%sortCriterion1 = sizeAssets;

%Considering the Book to Market
sortCriterion1 = BM;

%Considering betas, we have to resize the data of betas
 load('Betas.mat');
% betas =  nan( numMonths , numAssets);
% betas(:,1) = monthlyReturns(:,1);
% betas( size(monthlyReturns,1) -  size(betaCoeff,1) +1:end,2:end) = betaCoeff ;
% sortCriterion1 = BetasNew(:,2:end);

%Criterion 2 is always Momentum
sortCriterion2 = skewnessMatrix;

%Define Number of Categories for each Criterion
numOfCategoriesCriterion1 = 3;
numOfCategoriesCriterion2 = 3;

%Initialization of Return Matrices
equallyBDReturns=nan(numMonths,2);
valueBDReturns=nan(numMonths,2);

%Initialization of Matrix of Returns
equallyDecileBD = nan( numMonths, numOfCategoriesCriterion1+1,numOfCategoriesCriterion2+1);
valueDecileBD = equallyDecileBD;


for iInvestDate=60:numMonths-1
    %At first we sort based on Criterion 1
    [sortedCriterion1,indexSortedCriterion1] = sort( sortCriterion1(iInvestDate,:),2);
    
    %Define how I will treat Nulls of two Criteria
    %And the number of "Assets" in each quantile
    numOfAssetsOfMonth = sum( ~isnan(sortCriterion1(iInvestDate,:)) & ...
        ~isnan(sortCriterion2(iInvestDate,:)) & ~isnan(monthlyReturns( iInvestDate+1, :) )) ;
   
    numDecileCriterion1 = round(numOfAssetsOfMonth / numOfCategoriesCriterion1);
    numDecileCriterion2 = round(numDecileCriterion1 / numOfCategoriesCriterion2);
   
    %Sorting according to Criterion 2
    %Indexing for each Decile
    stepCriterion1 = numOfAssetsOfMonth;
    for j=numOfCategoriesCriterion1:-1:1
     if (j~=1)
        [sortedCriterion1,indexSortedCriterion2] = sort( sortCriterion2(iInvestDate,...
           indexSortedCriterion1( stepCriterion1 - numDecileCriterion1+1: stepCriterion1)),2);
        
            stepCriterion1 = stepCriterion1 - numDecileCriterion1;
        else
           
          [sortedCriterion1,indexSortedCriterion2] = sort( sortCriterion2(iInvestDate,indexSortedCriterion1( 1:numDecileCriterion1  )  ),2);
        
        end
    
        stepCriterion2 = length(indexSortedCriterion2) ;  
        for k =numOfCategoriesCriterion2:-1:1
        
         if (k~=1)
                         
           equallyDecileBD(iInvestDate+1,k,j) = nanmean(monthlyReturns(iInvestDate+1, ...
                indexSortedCriterion1(indexSortedCriterion2...
                 (stepCriterion2 - numDecileCriterion2+1: stepCriterion2))));
                 
          weight = sizeAssets(iInvestDate, indexSortedCriterion1( indexSortedCriterion2...
             ( stepCriterion2 - numDecileCriterion2+1: stepCriterion2  ) ))...
               /nansum(sizeAssets(iInvestDate, indexSortedCriterion1( indexSortedCriterion2...
                    ( stepCriterion2 - numDecileCriterion2+1: stepCriterion2  )  ) )); 
                
                valueDecileBD(iInvestDate+1,k,j) = nansum( monthlyReturns(iInvestDate+1, indexSortedCriterion1...
                    ( indexSortedCriterion2( stepCriterion2 - numDecileCriterion2+1: stepCriterion2 )  ) )  .*weight );

                stepCriterion2 = stepCriterion2 - numDecileCriterion2;
                
         else
             
     equallyDecileBD(iInvestDate+1,k,j) = nanmean(monthlyReturns(iInvestDate+1,  indexSortedCriterion1( indexSortedCriterion2( 1:numDecileCriterion2 ))));
               
        weight = sizeAssets(iInvestDate,  indexSortedCriterion1( indexSortedCriterion2...
           ( 1:numDecileCriterion2 ) ))...
            /nansum(sizeAssets(iInvestDate,  indexSortedCriterion1( indexSortedCriterion2...
           ( 1:numDecileCriterion2  )))); 

    valueDecileBD(iInvestDate+1,k,j) = nansum( monthlyReturns(iInvestDate+1, indexSortedCriterion1...
                 ( indexSortedCriterion2( 1:numDecileCriterion2 )  ) )  .*weight );
       
         end
                  
        end
             
        treatNan1 = sum(isnan( sortedCriterion1));
        position3= indexSortedCriterion1( indexSortedCriterion2...
              (length(indexSortedCriterion2) - treatNan1 - numDecileCriterion2+1: length(indexSortedCriterion2) - treatNan1 ));
      
        position1= indexSortedCriterion1( indexSortedCriterion2...
               (1:numDecileCriterion2));
        %Calculation for Equally weighted
        equallyDecileBD(iInvestDate+1,numOfCategoriesCriterion2+1,j)...
            =nanmean(monthlyReturns(iInvestDate+1,position3))-...
             nanmean(monthlyReturns(iInvestDate+1,position1));  

        %Calculation for Value weighted
 weightLong= sizeAssets(iInvestDate,position3)/sum(sizeAssets(iInvestDate,position3));
 weightShort= sizeAssets(iInvestDate,position1)/sum(sizeAssets(iInvestDate,position1));

    
    valueDecileBD(iInvestDate+1,numOfCategoriesCriterion2+1,j)...
   =sum(monthlyReturns(iInvestDate+1,position3).* weightLong)-...
    sum(monthlyReturns(iInvestDate+1,position1).*weightShort);      
                            
    end    
    
    for m=1:numOfCategoriesCriterion2+1
        equallyDecileBD(iInvestDate+1,m,numOfCategoriesCriterion1+1) = nanmean( squeeze(equallyDecileBD(iInvestDate+1,m,1:numOfCategoriesCriterion2)) ,1 );
        valueDecileBD(iInvestDate+1,m,numOfCategoriesCriterion1+1) = nanmean( squeeze(valueDecileBD(iInvestDate+1,m,1:numOfCategoriesCriterion2)) ,1 );
    end
        
    
end
%Ploting of Skewness Portfolio 
ewRwithNoNanBD=equallyDecileBD(:,numOfCategoriesCriterion2+1,numOfCategoriesCriterion1);ewRwithNoNanBD(isnan(ewRwithNoNanBD))=0;
vwRwithNoNanBD=valueDecileBD(:,numOfCategoriesCriterion2+1,numOfCategoriesCriterion1);vwRwithNoNanBD(isnan(vwRwithNoNanBD))=0;
 

dm=load('DATESIZE')  
figure
plot(datenum(num2str(dm.Date2(60:1098)),'yyyymm'),cumprod(1+ewRwithNoNanBD(60:1098)));
hold on
plot(datenum(num2str(dm.Date2(60:1098)),'yyyymm'),cumprod(1+vwRwithNoNanBD(60:1098)));
datetick

theAverageEqReturnBD = squeeze(nanmean( equallyDecileBD,1 ));
theAverageEqReturnBD = theAverageEqReturnBD*100;

theAverageVlReturnBD = squeeze(nanmean( valueDecileBD,1 ));
theAverageVlReturnBD = theAverageVlReturnBD*100;
% %Printing results in excel for size and skewness
%xlswrite('averageEwRSizeSkewness',theAverageEqReturnBD);
%xlswrite('averageVwRSizeSkewness',theAverageVlReturnBD);

% %Printing results in excel for book to market and skewness
 xlswrite('averageEwRBMSkew',theAverageEqReturnBD);
 xlswrite('averageVwRBMSkew',theAverageVlReturnBD);

%Initialization of Matrices
capmA = nan(1,numOfCategoriesCriterion1+1);
capmTstatisticsAlpha = capmA;
famaA =capmA;
famaTstatisticsAlpha = famaA;
capmTstatisticsAlphaV= famaA;
capmAV= famaA;
famaTstatisticsAlphaV= famaA;
famaAV= famaA;

onesTstatistics = nan(2,numOfCategoriesCriterion1+1);
onesRegression = ones( numMonths,1 );

% 
for j= 1:numOfCategoriesCriterion1+1
        
       %Treating Null in each iteration for appropriate regression
       treatNull = (isnan(equallyDecileBD(:,numOfCategoriesCriterion2+1,j)));
       
       %Definition of independent variables for FAMA FRENCH Regression
       famaIndependentVar = [MKT(~treatNull), SMB(~treatNull), HML(~treatNull)  ];
       
       %Regression with one's  
       [~,onesNwErrorsE, onesCoeffE ] =  hac( onesRegression(~treatNull),  equallyDecileBD(~treatNull,numOfCategoriesCriterion2+1,j)*100, 'intercept', false );
       onesTstatistics(1,j) = onesCoeffE(1)/onesNwErrorsE(1);
       [~,onesNwErrorsV, onesCoeffV ] =  hac( onesRegression(~treatNull),  valueDecileBD(~treatNull,numOfCategoriesCriterion2+1,j)*100, 'intercept', false );
       onesTstatistics(2,j) = onesCoeffV(1)/onesNwErrorsV(1);
       
       
       %CAPM Regression for Equally Weighted Portfolio
       %[covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF(~treatNull),  ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull) );
       [~,capmNwErrorsE, capmCoeffE ] =  hac( MKT(~treatNull),  equallyDecileBD(~treatNull,numOfCategoriesCriterion2+1,j)*100 );
       capmTstatisticsAlpha(j) = capmCoeffE(1)/capmNwErrorsE(1);
       capmA(j) = capmCoeffE(1);
       
       %Fama–MacBeth Regression for Equally Weighted Portfolio
       %[covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull)  );
       [~,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, equallyDecileBD(~treatNull,numOfCategoriesCriterion2+1,j)*100   );
       famaTstatisticsAlpha(j) = famaCoeffE(1)/famaNwErrorsE(1);
       famaA(j) = famaCoeffE(1);
       
       %CAPM Regression for Value Weighted Portfolio
       %[covarCapmV,capmNwErrorsV, capmCoeffV ] =  hac( marketMinusRF(~treatNull),  vwROfDecile(~treatNull,6,j)*100 - RF(~treatNull) );
       [~,capmNwErrorsV, capmCoeffV ] =  hac( MKT(~treatNull),  valueDecileBD(~treatNull,numOfCategoriesCriterion2+1,j)*100 );
       capmTstatisticsAlphaV(j) = capmCoeffV(1)/capmNwErrorsV(1);
       capmAV(j) = capmCoeffV(1);
       
       %Fama–MacBeth Regression for  Value Weighted Portfolio
       %[covarFAMAV,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar, vwROfDecile(~treatNull,6,j)*100 - RF(~treatNull)  );
       [~,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar, valueDecileBD(~treatNull,numOfCategoriesCriterion2+1,j)*100   );
       famaTstatisticsAlphaV(j) = famaCoeffV(1)/famaNwErrorsV(1);
       famaAV(j) = famaCoeffV(1);
       
       
end

%Printing results for size and skewness depending sorting
%resultsforExcelequally = [theAverageEqReturnBD ; onesTstatistics(1,:) ; capmA ; capmTstatisticsAlpha ; famaA ; famaTstatisticsAlpha;];
 % resultsforExcelvalue = [theAverageVlReturnBD ; onesTstatistics(2,:) ; capmAV ; capmTstatisticsAlphaV ; famaAV ; famaTstatisticsAlphaV;   ];
%xlswrite('BivariateDependSortEquallySizeSkew',resultsforExcelequally);
%xlswrite('BivariateDependSortValueSizeSkew',resultsforExcelvalue);

% %Printing results for book to market and skewness depending sorting
 resultsforExcelequallyBMSkew = [theAverageEqReturnBD ; onesTstatistics(1,:) ; capmA ; capmTstatisticsAlpha ; famaA ; famaTstatisticsAlpha;];
 resultsforExcelvalueBMSkew = [theAverageVlReturnBD ; onesTstatistics(2,:) ; capmAV ; capmTstatisticsAlphaV ; famaAV ; famaTstatisticsAlphaV;   ];
 xlswrite('BivariateDependSortEquallyBMSkew',resultsforExcelequallyBMSkew);
 xlswrite('BivariateDependSortValueBMSkew',resultsforExcelvalueBMSkew );
 


             