%Defining Sorting Criteria

%For Size
%sortCriterion1 = sizeAssets;

%For BM
%load('BM.mat');
%BM = [ monthlyReturns(:,1) ,BM(1:1098,:)  ];
%sortCriterion1 = BM;

% For Betas
%  load('Betas.mat');
%  BetasNew =  nan( numMonths , numAssets+1);
%  BetasNew(:,1) = monthlyReturns(:,1);
%  BetasNew( size(monthlyReturns,1) -  size(betaCoeff,1)+1:end,2:end) = betaCoeff ;
    sortCriterion1 = BetasNew(:,2:end);
 
%Criterion 2 is always Momentum
sortCriterion2 = skewnessMatrix;

%Define Number of Categories for each Criterion
numOfCategoriesCriterion1 = 3;
numOfCategoriesCriterion2 = 3;


%Main Part for the Estimation of Return
%For InDependent Sorts


%Initialization of Return Matrices
ewRIND=nan(numMonths,4);
vwRIND=nan(numMonths,4);


%Initialization of Matrix of Returns of 9*9 matrix
ewROfQuintile = nan( numMonths, numOfCategoriesCriterion1+1);
vwROfQuintile = ewROfQuintile;

for i=60:numMonths-1
    
    numOfAssetsOfMonth = sum( ~isnan(sortCriterion1(i,:)) & ...
        ~isnan(sortCriterion2(i,:)) & ~isnan(monthlyReturns( i+1,:) )) ;  
     
   numAssetsInQuantileCriterion1 = round(numOfAssetsOfMonth / numOfCategoriesCriterion1);
   numAssetsInQuantileCriterion2 = round(numOfAssetsOfMonth / numOfCategoriesCriterion2);

   
    %Sorting according to Criterion 1
    [sortedMatrixofReturnsCriterion1,indexSortedCriterion1] = sort( sortCriterion1(i,:),2);
    [sortedMatrixofReturnsCriterion2,indexSortedCriterion2] = sort( sortCriterion2(i,:),2);
       
    %Sorting according to Criterion 2
    %Indexing for each Decile
    stepCriterion1 = numOfAssetsOfMonth;
    
    for j1 = numOfCategoriesCriterion1:-1:1
            
        if (j1~=1)
            
            positionInIndex1 = stepCriterion1 - numAssetsInQuantileCriterion1+1: stepCriterion1;
            
        else

            positionInIndex1 =  1:numAssetsInQuantileCriterion1;
            
        end
        stepCriterion2 = numOfAssetsOfMonth;
        for j2 = numOfCategoriesCriterion2:-1:1
            %Estimation of Returns for Equally Weighted Portfolio
              
            if(j2~=1)
                
            positionInIndex2 = stepCriterion2 - numAssetsInQuantileCriterion2+1: stepCriterion2;
                
            else
                
            positionInIndex2 =  1:numAssetsInQuantileCriterion2;
            
            end
            
            indexAssetsOfCategory = intersect( indexSortedCriterion1( positionInIndex1 ) , indexSortedCriterion2(positionInIndex2) );
            
            ewROfQuintile(j1,j2) = nanmean(monthlyReturns(i+1, ...
                     indexAssetsOfCategory));
            
            %Estimation of Returns for Value Weighted Portfolio
%             weight = sizeAssets(i, 1 + indexAssetsOfCategory( positionInIndex ))...
%                     /nansum(sizeAssets(i, 1+ indexAssetsOfCategory( positionInIndex ))) ; 
% 
%             vwROfDecile(i+1,k,j) = nansum( monthlyReturns(i+1, 1 +  indexAssetsOfCategory( positionInIndex )).*weight );
%                 
            
            weight = sizeAssets(i, indexAssetsOfCategory)...
                    /nansum(sizeAssets(i, indexAssetsOfCategory )) ; 

            vwROfQuintile(j1,j2) = nansum( monthlyReturns(i+1,  ...
                indexAssetsOfCategory).*weight );
                
            stepCriterion2 = stepCriterion2 - numAssetsInQuantileCriterion2;
        

        end
        
        stepCriterion1 = stepCriterion1 - numAssetsInQuantileCriterion1;
    
         ewRIND(i,j1) = ewROfQuintile(j1,3) - ewROfQuintile(j1,1);
         vwRIND(i,j1) = vwROfQuintile(j1,3) - vwROfQuintile(j1,1);

        
    end
        
    ewRIND(i,4) = nanmean(ewRIND(i,1:3));
    vwRIND(i,4) = nanmean(vwRIND(i,1:3));
   
    
end

averageEwRIND = squeeze(nanmean( ewRIND,1 ));
averageEwRIND = averageEwRIND*100;


averageVwRIND = squeeze(nanmean( vwRIND,1 )); 
averageVwRIND = averageVwRIND*100;
% 
% 
% %% Write in Excel the Results
% xlswrite('averageEwRSizeMomentum',averageEwR);
% xlswrite('averageVwRSizeMomentum',averageVwR);
% 
% % xlswrite('averageEwR_BM_Momentum',averageEwR);
% % xlswrite('averageVwR_BM_Momentum',averageVwR);
% 
% % xlswrite('averageEwR_Betas_Momentum',averageEwR);
% % xlswrite('averageVwR_Betas_Momentum',averageVwR);
% %%
% 
% %%% Regression Part %%%%

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
       treatNull = (isnan(ewRIND(:,j)));
       
       %Definition of independent variables for FAMA FRENCH Regression
       famaIndependentVar = [ MKT(~treatNull), SMB(~treatNull), HML(~treatNull)  ];
       
       
       %Regression with one's  
       [~,onesNwErrorsE, onesCoeffE ] =  hac( onesRegression(~treatNull),  ewRIND(~treatNull,j)*100, 'intercept', false );
       onesTstatistics(1,j) = onesCoeffE(1)/onesNwErrorsE(1);
       [~,onesNwErrorsV, onesCoeffV ] =  hac( onesRegression(~treatNull),  vwRIND(~treatNull,j)*100, 'intercept', false );
       onesTstatistics(2,j) = onesCoeffV(1)/onesNwErrorsV(1);
       
       
       %CAPM Regression for Equally Weighted Portfolio
       %[covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF(~treatNull),  ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull) );
       [~,capmNwErrorsE, capmCoeffE ] =  hac( MKT(~treatNull),  ewRIND(~treatNull,j)*100 );
       capmTstatisticsAlpha(j) = capmCoeffE(1)/capmNwErrorsE(1);
       capmA(j) = capmCoeffE(1);
       
       %Fama–MacBeth Regression for Equally Weighted Portfolio
       %[covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull)  );
       [~,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, ewRIND(~treatNull,j)*100   );
       famaTstatisticsAlpha(j) = famaCoeffE(1)/famaNwErrorsE(1);
       famaA(j) = famaCoeffE(1);
       
       %CAPM Regression for Value Weighted Portfolio
       
       [~,capmNwErrorsV, capmCoeffV ] =  hac( MKT(~treatNull),  vwRIND(~treatNull,j)*100 );
       capmTstatisticsAlphaV(j) = capmCoeffV(1)/capmNwErrorsV(1);
       capmAV(j) = capmCoeffV(1);
       
       %Fama–MacBeth Regression for  Value Weighted Portfolio
       [~,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar, vwRIND(~treatNull,j)*100   );
       famaTstatisticsAlphaV(j) = famaCoeffV(1)/famaNwErrorsV(1);
       famaAV(j) = famaCoeffV(1);
       
       
end


% equallySizeSkew = [averageEwRIND ; onesTstatistics(1,:) ; capmA ; capmTstatisticsAlpha ; famaA ; famaTstatisticsAlpha; ];
% xlswrite('BivariateIndependentSortingEquallySizeSkew',equallySizeSkew);
% 
% valueSizeSkewInd = [averageVwRIND ; onesTstatistics(2,:) ; capmAV ; capmTstatisticsAlphaV ; famaAV ; famaTstatisticsAlphaV; ];
% xlswrite('BivariateIndependentSortingValueSizeSkew',valueSizeSkewInd);

% equallyBMSkew = [averageEwRIND ; onesTstatistics(1,:) ; capmA ; capmTstatisticsAlpha ; famaA ; famaTstatisticsAlpha; ];
% xlswrite('BivariateIndependentSortingEquallyBMSkew',equallyBMSkew);
% 
% valueBMSkewInd = [averageVwRIND ; onesTstatistics(2,:) ; capmAV ; capmTstatisticsAlphaV ; famaAV ; famaTstatisticsAlphaV; ];
% xlswrite('BivariateIndependentSortingValueBMSkew',valueBMSkewInd);

equallyBetasSkew = [averageEwRIND ; onesTstatistics(1,:) ; capmA ; capmTstatisticsAlpha ; famaA ; famaTstatisticsAlpha; ];
xlswrite('BivariateIndependentSortingEquallyBetasSkew',equallyBetasSkew);

valueBetasSkewInd = [averageVwRIND ; onesTstatistics(2,:) ; capmAV ; capmTstatisticsAlphaV ; famaAV ; famaTstatisticsAlphaV; ];
xlswrite('BivariateIndependentSortingValueBetasSkew',valueBetasSkewInd);