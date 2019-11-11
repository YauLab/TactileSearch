function ompRes = MinBICSplitTPTA(ompRes)
% Modified version of selecting models based on BIC - now finding BIC for
% target present and target absent trials separately.

nRes = length(ompRes);
bPlot = false;
for ii = 1:nRes

    RT = ompRes(ii).res.response.RT;
    bTP = ompRes(ii).res.params.iTarget > 0;
    
    [D,~] = Make2ndPBSDict(ompRes(ii).res); D = double(D);
    nTrials = length(RT);
    nTerms = size(D,2);
    
    %% 17-term fit (linear terms only):
    w17 = FindNWeights(ompRes(ii).res,1:17,D);
    w17n = zeros([1 nTerms]);
    w17n(1:17) = w17;
    [RMSE_17,R2_17,yhat_17] = FitNLinReg3(ompRes(ii).res,w17n,D,bPlot);
    ompRes(ii).linear.R2 = R2_17;
    ompRes(ii).linear.terms = 1:17;
    ompRes(ii).linear.coefs = w17n;
    ompRes(ii).linear.RMSE = RMSE_17;
    ompRes(ii).linear.RMSE = RMSE_17;
    
    % Find loglikelihood of this model:
    muHatTP = mean(yhat_17(bTP)); muHatTA = mean(yhat_17(~bTP));
    sigmaHatTP = std(yhat_17(bTP)); sigmaHatTA = std(yhat_17(~bTP));
    NLLtp = normlike([muHatTP sigmaHatTP],RT(bTP));
    NLLta = normlike([muHatTA sigmaHatTA],RT(~bTP));
    negloglik17 = NLLtp + NLLta;
    ompRes(ii).linear.NLL = negloglik17;
    
    % Find AIC & BIC for this linear fit:
    nLinTerms = 17;
    ompRes(ii).linear.AIC = 2*negloglik17 + 2*nLinTerms;
    ompRes(ii).linear.BIC = 2*negloglik17 + log(nTrials)*nLinTerms;
    
    R2 = nan([1 nTerms]); RMSE = nan([1 nTerms]);
    NLL = nan([1 nTerms]); NLLtp = nan([1 nTerms]); NLLta = nan([1 nTerms]);
    AIC = nan([1 nTerms]); BIC = nan([1 nTerms]);
    
    for nT1 = 17:nTerms
        %% Find NLL & R2 for each model:
        
        % Find coefs for this set of terms:
        w = ompRes(ii).coefs(nT1,:);
        
        % Fit model to find yhat:
        [RMSE(nT1),R2(nT1),yhat] = FitNLinReg3(ompRes(ii).res,w,D,bPlot);
        
        % NLL    (separate TP & TA distributions)
        muHatTP = mean(yhat(bTP)); muHatTA = mean(yhat(~bTP));
        sigmaHatTP = std(yhat(bTP)); sigmaHatTA = std(yhat(~bTP));
        NLLtp(nT1) = normlike([muHatTP sigmaHatTP],RT(bTP));
        NLLta(nT1) = normlike([muHatTA sigmaHatTA],RT(~bTP));
        NLL(nT1) = NLLtp(nT1) + NLLta(nT1);
        
        % Find AIC & BIC:
        AIC(nT1) = 2*NLL(nT1) + 2*nT1;
        BIC(nT1) = 2*NLL(nT1) + log(nTrials)*nT1;
    end
%     figure('pos',[0 0 1200 600]);yyaxis left;plot(BIC,'b');yyaxis right;plot(RMSE,'r');
    %% Find best models & store results:
    % Based on optimizing NLL:
    [minNLL,iNLL] = nanmin(NLL);
    terms1 = ompRes(ii).terms{iNLL};
    ompRes(ii).optNLL.nTerms = length(terms1);
    ompRes(ii).optNLL.nLinTerms = sum(terms1>1 & terms1<=17);
    ompRes(ii).optNLL.terms = terms1;
    ompRes(ii).optNLL.coefs = ompRes(ii).coefs(iNLL,:);
    ompRes(ii).optNLL.R2 = R2(iNLL);
    ompRes(ii).optNLL.RMSE = RMSE(iNLL);
    ompRes(ii).optNLL.NLL = minNLL;
    ompRes(ii).optNLL.AIC = AIC(iNLL);
    ompRes(ii).optNLL.BIC = BIC(iNLL);
    % Based on optimizing AIC:
    [minAIC,iAIC] = nanmin(AIC);
    terms1 = ompRes(ii).terms{iAIC};
    ompRes(ii).optAIC.nTerms = length(terms1);
    ompRes(ii).optAIC.nLinTerms = sum(terms1>1 & terms1<=17);
    ompRes(ii).optAIC.terms = terms1;
    ompRes(ii).optAIC.coefs = ompRes(ii).coefs(iAIC,:);
    ompRes(ii).optAIC.R2 = R2(iAIC);
    ompRes(ii).optAIC.RMSE = RMSE(iAIC);
    ompRes(ii).optAIC.NLL = NLL(iAIC);
    ompRes(ii).optAIC.AIC = minAIC;
    ompRes(ii).optAIC.BIC = BIC(iAIC);
    % Based on optimizing BIC:
    [minBIC,iBIC] = nanmin(BIC);
    terms1 = ompRes(ii).terms{iBIC};
    ompRes(ii).optBIC.nTerms = length(terms1);
    ompRes(ii).optBIC.nLinTerms = sum(terms1>1 & terms1<=17);
    ompRes(ii).optBIC.terms = terms1;
    ompRes(ii).optBIC.coefs = ompRes(ii).coefs(iBIC,:);
    ompRes(ii).optBIC.R2 = R2(iBIC);
    ompRes(ii).optBIC.RMSE = RMSE(iBIC);
    ompRes(ii).optBIC.NLL = NLL(iBIC);
    ompRes(ii).optBIC.AIC = AIC(iBIC);
    ompRes(ii).optBIC.BIC = minBIC;
    
    %% Store general results:
    ompRes(ii).R2 = R2;
    ompRes(ii).RMSE = RMSE;
    ompRes(ii).NLL = NLL;
    ompRes(ii).AIC = AIC;
    ompRes(ii).BIC = BIC;
end

end