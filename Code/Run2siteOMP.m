function results = Run2siteOMP(resAll,subjIDs,bFixLin,bSave)
% results = Run2siteOMP(resAll,subjIDs,bFixLin,bSave)
%
% Run OMP for all (137) terms, for nIter1 number of
% iterations. Store NLL for each iteration & find which sequence occurs
% most frequently across the iterations.
%
% ejh 04/08/19

%% Initialize:
if ~exist('bSave','var'), bSave = true; end
if ~exist('bFixLin','var'), bFixLin = true; end
if ~exist('subjIDs','var'), subjIDs = 1:length(resAll); end
nIter = 10000;                   % 100it = 20-30s, 1000it = 3-4m
nSubj = length(subjIDs);
ex1 = struct('res',[],'D',[],'terms',[],'coefs',[],'R2',[],'NLL',[],'AIC',[],'BIC',[]);
ex2 = struct('nTerms',[],'nLinTerms',[],'terms',[],'coefs',[],'R2',[],'NLL',[],'AIC',[],'BIC',[]);
results(max(subjIDs),1) = ex1;
for ii = 1:max(subjIDs)
    results(ii).linear = ex2;
    results(ii).optNLL = ex2;
    results(ii).optAIC = ex2;
    results(ii).optBIC = ex2;
end

if nSubj > 1
    tFreq = num2str(resAll{subjIDs(1)}.params.tFreq);
    if any(cellfun(@(x) any(x>0),resAll{subjIDs(1)}.params.dOffset))
        psStr = 'PhSc';
    else
        psStr = 'Sync';
    end
else
    tFreq = num2str(resAll(subjIDs(1)).params.tFreq);
    if any(cellfun(@(x) any(x>0),resAll(subjIDs(1)).params.dOffset))
        psStr = 'PhSc';
    else
        psStr = 'Sync';
    end
end

if bFixLin
    keepTerms = 1:17;
else
    keepTerms = [];
end

%%
for pp = 1:nSubj
    iSubj = subjIDs(pp);
    if nSubj > 1
        res = resAll{iSubj};
    else
        res = resAll(iSubj);
    end
    
    results(iSubj).res = res;
    RT = res.response.RT;
    bTP = res.params.iTarget > 0;
    [D,NLsites] = Make2ndPBSDict(res); D = double(D);
    results(iSubj).D = D;
    nTrials = length(RT);
    nTerms = size(D,2);
    
    %% 17-term fit (linear terms only):
    bPlot = false;
    w17 = FindNWeights(res,1:17,D);
    w17n = zeros([1 nTerms]);
    w17n(1:17) = w17;
    [RMSE_17,R2_17,simRT_17] = FitNLinReg3(res,w17n,D,bPlot); 
    results(iSubj).linear.R2 = R2_17;
    results(iSubj).linear.terms = 1:17;
    results(iSubj).linear.coefs = w17n;
    
    % Find loglikelihood of this model:
    sqErr17 = (RT'-simRT_17).^2;
    pd17 = fitdist(sqErr17,'Normal');
    
    % Find likelihood:
    negloglik17 = negloglik(pd17);
    
    results(iSubj).linear.NLL = negloglik17;
    
    % Find AIC & BIC for this linear fit:
    results(iSubj).linear.AIC = 2*negloglik17 + 2*17;
    results(iSubj).linear.BIC = 2*negloglik17 + log(nTrials)*17;
    
    %% Run OMP (nIter1 iterations):
    [x1,b1] = OMP_ThakurMod(D,RT',nTerms,keepTerms,nIter);
        
    %% For each number of terms, find the coefficients:
    coefs = nan([nTerms nTerms]); % row-wise
    mTerms = cell([1 nTerms]);
    R2 = nan([1 nTerms]); NLL = nan([1 nTerms]);
    AIC = nan([1 nTerms]); BIC = nan([1 nTerms]);

    % Fill in linear coefficients:
    blank = zeros([1 nTerms]);
    for nT = 1:nTerms
        % Find most likely set of nTerms terms:
        b = b1(1:nT,:);
        ct1 = tabulate(b(:));

        [~,ts] = maxk(ct1(:,2),nT);
        mTerms{nT} = ts;
        
        % Find coefficients & store in full coefs array:
        % --> what to do with constant if not first term??
        D1 = D(:,ts);
        weights = D1 \ RT'; % based on residual signal
        w = blank;
        w(ts) = weights;
        coefs(nT,:) = w;
        
        %% Find NLL & R2 for each model:
        [~,R2(nT),yhat] = FitNLinReg3(res,w,D,bPlot);
        
        % NLL    (separate TP & TA distributions)
        muHatTP = mean(yhat(bTP)); muHatTA = mean(yhat(~bTP));
        sigmaHatTP = std(yhat(bTP)); sigmaHatTA = std(yhat(~bTP));
        
        NLLtp = normlike([muHatTP sigmaHatTP],RT(bTP));
        NLLta = normlike([muHatTA sigmaHatTA],RT(~bTP));
        NLL(nT1) = NLLtp + NLLta;

        % Find AIC & BIC:
        AIC(nT) = 2*NLL(nT) + 2*nT;
        BIC(nT) = 2*NLL(nT) + log(nTrials)*nT;
    end
    
    %% Find best models & store results:
    % Based on optimizing NLL:
    [minNLL,iNLL] = nanmin(NLL);
    terms1 = mTerms{iNLL};
    results(iSubj).optNLL.nTerms = length(terms1);
    results(iSubj).optNLL.nLinTerms = sum(terms1>1 & terms1<=17);
    results(iSubj).optNLL.terms = terms1;
    results(iSubj).optNLL.coefs = coefs(iNLL,:);
    results(iSubj).optNLL.R2 = R2(iNLL);
    results(iSubj).optNLL.NLL = minNLL;
    results(iSubj).optNLL.AIC = AIC(iNLL);
    results(iSubj).optNLL.BIC = BIC(iNLL);
    % Based on optimizing AIC:
    [minAIC,iAIC] = nanmin(AIC);
    terms1 = mTerms{iAIC};
    results(iSubj).optAIC.nTerms = length(terms1);
    results(iSubj).optAIC.nLinTerms = sum(terms1>1 & terms1<=17);
    results(iSubj).optAIC.terms = terms1;
    results(iSubj).optAIC.coefs = coefs(iAIC,:);
    results(iSubj).optAIC.R2 = R2(iAIC);
    results(iSubj).optAIC.NLL = NLL(iAIC);
    results(iSubj).optAIC.AIC = minAIC;
    results(iSubj).optAIC.BIC = BIC(iAIC);
    % Based on optimizing BIC:
    [minBIC,iBIC] = nanmin(BIC);
    terms1 = mTerms{iBIC};
    results(iSubj).optBIC.nTerms = length(terms1);
    results(iSubj).optBIC.nLinTerms = sum(terms1>1 & terms1<=17);
    results(iSubj).optBIC.terms = terms1;
    results(iSubj).optBIC.coefs = coefs(iBIC,:);
    results(iSubj).optBIC.R2 = R2(iBIC);
    results(iSubj).optBIC.NLL = NLL(iBIC);
    results(iSubj).optBIC.AIC = AIC(iBIC);
    results(iSubj).optBIC.BIC = minBIC;
    
    %% Store general results:
    results(iSubj).terms = mTerms;
    results(iSubj).coefs = coefs;
    results(iSubj).R2 = R2;
    results(iSubj).NLL = NLL;
    results(iSubj).AIC = AIC;
    results(iSubj).BIC = BIC;
    
    %% Plot results:
    mSize = 10; lWidth = 1.5;        
    namestr = [res.params.initials ', ' tFreq 'T, ' psStr];
    figure('color',[1 1 1],'pos',[0 0 1000 1000],'Name',namestr);
    x = 1:nTerms;
    cmap = cbrewer('qual','Set1',3);
    
    % Plot NLL:
    subplot(2,2,1); hold on;
    plot(x,NLL,'color',cmap(1,:),'LineWidth',2); box off;
    plot(iNLL,results(iSubj).optNLL.NLL,'ko','MarkerSize',mSize,'LineWidth',lWidth);
    plot(17,negloglik17,'r*','MarkerSize',mSize);
    xlabel('# of Terms'); ylabel('Negative Log Likelihood');
    set(gca,'FontSize',18);
    title(['NLL (' num2str(iNLL) ' Terms, ' num2str(results(iSubj).optNLL.nLinTerms) ' Linear)']);
    roundY = 50;
    minY = floor(nanmin([NLL negloglik17])/roundY)*roundY;
    maxY = ceil(nanmax([NLL negloglik17])/roundY)*roundY;
    axis([0 nTerms minY maxY]);
    
    % Plot AIC:
    subplot(2,2,2); hold on;
    plot(x,AIC,'color',cmap(2,:),'LineWidth',2); box off;
    plot(iAIC,results(iSubj).optAIC.AIC,'ko','MarkerSize',mSize,'LineWidth',lWidth);
    plot(17,results(iSubj).linear.AIC,'r*','MarkerSize',mSize);
    xlabel('# of Terms'); ylabel('AIC');
    set(gca,'FontSize',18);
    title(['AIC (' num2str(iAIC) ' Terms, ' num2str(results(iSubj).optAIC.nLinTerms) ' Linear)']);
    roundY = 100;
    minY = floor(nanmin([AIC results(iSubj).linear.AIC])/roundY)*roundY;
    maxY = ceil(nanmax([AIC results(iSubj).linear.AIC])/roundY)*roundY;
    axis([0 nTerms minY maxY]);
    
    % Plot BIC:
    subplot(2,2,3); hold on;
    plot(x,BIC,'color',cmap(3,:),'LineWidth',2); box off;
    plot(iBIC,results(iSubj).optBIC.BIC,'ko','MarkerSize',mSize,'LineWidth',lWidth);
    plot(17,results(iSubj).linear.BIC,'r*','MarkerSize',mSize);
    xlabel('# of Terms'); ylabel('BIC');
    set(gca,'FontSize',18);
    title(['BIC (' num2str(iBIC) ' Terms, ' num2str(results(iSubj).optBIC.nLinTerms) ' Linear)']);
    minY = floor(nanmin([BIC results(iSubj).linear.BIC])/roundY)*roundY;
    maxY = ceil(nanmax([BIC results(iSubj).linear.BIC])/roundY)*roundY;
    axis([0 nTerms minY maxY]);
    
    % Plot R2:
    subplot(2,2,4); hold on;
    plot(x,R2,'k'); box off;
    plot(x,repmat(R2(iNLL),size(x)),':','color',cmap(1,:),'LineWidth',2);
    plot(x,repmat(R2(iAIC),size(x)),'-.','color',cmap(2,:),'LineWidth',2);
    plot(x,repmat(R2(iBIC),size(x)),'--','color',cmap(3,:),'LineWidth',2);
    plot(17,R2_17,'r*');
    xlabel('# of Terms'); ylabel('R^{2}');
    set(gca,'FontSize',18);
    title('R^{2}');
    
    %% Save results:
    if bSave
        if bFixLin
            saveName  = ['sAll_' tFreq 'T_' psStr '_FixLinOMP_10000iter'];
        else
            saveName  = ['sAll_' tFreq 'T_' psStr '_FullLinOMP_10000iter'];
        end
        save(saveName,'results');
    end
end
end