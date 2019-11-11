function results = FitRTvSS(res,bToss)
% results = FitRTvSS(res,bToss)
%
% Find linear & quadratic fits plus their log likelihoods.
%
% ejh 03/24/19

%% Initialize:
if ~exist('bToss','var'), bToss = true; end
results = struct();
if ~bToss
    bCorrect = true(size(res.response.correct'));
else
    bCorrect = res.response.correct';    % correct responses
end
bTarget = res.params.iTarget > 0;    % target-present trials

% Change ND to set size:
setSize = res.params.ND;
setSize(bTarget) = setSize(bTarget) + 1;

nSetSize = 1:8;
nSS = length(nSetSize);

%% Collect data:
% Pre-allocate:
RT_CD = cell([1 nSS]);               % store all RTs, split by SS
RT_CR = cell([1 nSS]);
mRT_CD = nan([1 nSS]);               % mean @ each SS
mRT_CR = nan([1 nSS]);
stdRT_CD = nan([1 nSS]);             % std dev @ each SS
stdRT_CR = nan([1 nSS]);
RT_CD1 = []; RT_CR1 = [];            % store all RTs, unsplit
x_CD = []; x_CR = [];
for jj = 1:nSS
    % CORRECT DETECT
    bSS_CD = setSize == jj & bCorrect & bTarget;
    RT = res.response.RT(bSS_CD);
    mRT_CD(jj) = nanmean(RT);        % mean RT
    stdRT_CD(jj) = nanstd(RT);       % std RT
    % Store in full set of RTs:
    RT_CD{jj} = RT;
    RT_CD1 = [RT_CD1 RT];
    x_CD = [x_CD repmat(jj,size(RT))];
    
    % CORRECT REJECT
    bSS_CR = setSize == jj & bCorrect & ~bTarget;
    RT = res.response.RT(bSS_CR);
    mRT_CR(jj) = nanmean(RT);
    stdRT_CR(jj) = nanstd(RT);
    % Store in full set of RTs:
    RT_CR{jj} = RT;
    RT_CR1 = [RT_CR1 RT];
    x_CR = [x_CR repmat(jj,size(RT))];
end

%% Linear fits:
x = nSetSize; y = mRT_CD; bNan = isnan(y);            % CD, point fit
x = x(~bNan); y = y(~bNan);
[p1,S1] = polyfit(x,y,1);
mLinCD = p1(1); bLinCD = p1(2);
data1 = polyval(p1,nSetSize,S1); % simulate data

x = x_CD; y = RT_CD1; bNan = isnan(y);                % CD, cloud fit
x = x(~bNan); y = y(~bNan);
[p1c,S1c] = polyfit(x,y,1);
mLinCDc = p1c(1); bLinCDc = p1c(2);
data1c = polyval(p1c,nSetSize,S1c); % simulate data

x = nSetSize; y = mRT_CR; bNan = isnan(y);            % CR, point fit
x = x(~bNan); y = y(~bNan);
[p2,S2] = polyfit(x,y,1);
mLinCR = p2(1); bLinCR = p2(2);
data2 = polyval(p2,nSetSize,S2); % simulate data

x = x_CR; y = RT_CR1; bNan = isnan(y);
x = x(~bNan); y = y(~bNan);
[p2c,S2c] = polyfit(x,y,1);                          % CR, cloud fit
mLinCRc = p2c(1); bLinCRc = p2c(2);
data2c = polyval(p2c,nSetSize,S2c); % simulate data

%% Quadratic fits:
x = nSetSize; y = mRT_CD; bNan = isnan(y);           % CD, point fit
x = x(~bNan); y = y(~bNan);
[p3,S3] = polyfit(x,y,2);
m2QuadCD = p3(1); mQuadCD = p3(2); bQuadCD = p3(3);
data3 = polyval(p3,nSetSize,S3); % simulate data

x = x_CD; y = RT_CD1; bNan = isnan(y);               % CD, cloud fit
x = x(~bNan); y = y(~bNan);
[p3c,S3c] = polyfit(x,y,2);
m2QuadCDc = p3c(1); mQuadCDc = p3c(2); bQuadCDc = p3c(3);
data3c = polyval(p3c,nSetSize,S3c); % simulate data

x = nSetSize; y = mRT_CR; bNan = isnan(y);           % CR, point fit
x = x(~bNan); y = y(~bNan);
[p4,S4] = polyfit(x,y,2);
m2QuadCR = p4(1); mQuadCR = p4(2); bQuadCR = p4(3);
data4 = polyval(p4,nSetSize,S4); % simulate data

x = x_CR; y = RT_CR1; bNan = isnan(y);               % CR, cloud fit
x = x(~bNan); y = y(~bNan);
[p4c,S4c] = polyfit(x,y,2);
m2QuadCRc = p4c(1); mQuadCRc = p4c(2); bQuadCRc = p4c(3);
data4c = polyval(p4c,nSetSize,S4c); % simulate data

%% Find residuals:
residCD=[]; residCDc=[]; residCR=[]; residCRc=[];
residqCD=[]; residqCDc=[]; residqCR=[]; residqCRc=[];
for ii = 1:nSS
    rCDlin = (RT_CD{ii} - data1(ii)).^2;         % linear, CD
    residCD = [residCD rCDlin];
    rCDc = (RT_CD{ii} - data1c(ii)).^2;       % linear cloud, CD
    residCDc = [residCDc rCDc];
    
    rCRlin = (RT_CR{ii} - data2(ii)).^2;         % linear, CR
    residCR = [residCR rCRlin];
    rCRc = (RT_CR{ii} - data2c(ii)).^2;       % linear cloud, CR
    residCRc = [residCRc rCRc];
    
    rCDquad = (RT_CD{ii} - data3(ii)).^2;         % quadratic, CD
    residqCD = [residqCD rCDquad];
    rCDc = (RT_CD{ii} - data3c(ii)).^2;       % quadratic cloud, CD
    residqCDc = [residqCDc rCDc];
    
    rCRquad = (RT_CR{ii} - data4(ii)).^2;         % quadratic, CR
    residqCR = [residqCR rCRquad];
    rCRc = (RT_CR{ii} - data4c(ii)).^2;       % quadratic cloud, CR
    residqCRc = [residqCRc rCRc]; %#ok<*AGROW>
end

%% Fit normal distribution to residuals to find variance in y:
if ~all(isnan(residCD)) && length(residCD)>2
    N_CD = fitdist(residCD','Normal');         % linear, CD
    NLL_CD = negloglik(N_CD);
else
    N_CD.mu=nan;N_CD.sigma=nan;
    NLL_CD = nan;
end
if ~all(isnan(residCDc)) && length(residCDc)>2
    N_CDc = fitdist(residCDc','Normal');    % linear cloud, CD
    NLL_CDc = negloglik(N_CDc);
else
    N_CDc.mu=nan;N_CDc.sigma=nan;
    NLL_CDc = nan;
end
if ~all(isnan(residCR)) && length(residCR)>2
    N_CR = fitdist(residCR','Normal');      % linear, CR
    NLL_CR = negloglik(N_CR);
else
    N_CR.mu=nan;N_CR.sigma=nan;
    NLL_CR = nan;
end
if ~all(isnan(residCRc)) && length(residCRc)>2
    N_CRc = fitdist(residCRc','Normal');    % linear cloud, CR
    NLL_CRc = negloglik(N_CRc);
else
    N_CRc.mu=nan;N_CRc.sigma=nan;
    NLL_CRc = nan;
end
if ~all(isnan(residqCD)) && length(residqCD)>2
    Nq_CD = fitdist(residqCD','Normal');    % quadratic, CD
    NLLq_CD = negloglik(Nq_CD);
else
    Nq_CD.mu=nan;Nq_CD.sigma=nan;
    NLLq_CD = nan;
end
if ~all(isnan(residqCDc)) && length(residqCDc)>2
    Nq_CDc = fitdist(residqCDc','Normal');  % quadratic cloud, CD
    NLLq_CDc = negloglik(Nq_CDc);
else
    Nq_CDc.mu=nan;Nq_CDc.sigma=nan;
    NLLq_CDc = nan;
end
if ~all(isnan(residqCR)) && length(residqCR)>2
    Nq_CR = fitdist(residqCR','Normal');    % quadratic, CR
    NLLq_CR = negloglik(Nq_CR);
else
    Nq_CR.mu=nan;Nq_CR.sigma=nan;
    NLLq_CR = nan;
end
if ~all(isnan(residqCRc)) && length(residqCRc)>2
    Nq_CRc = fitdist(residqCRc','Normal');  % quadratic cloud, CR
    NLLq_CRc = negloglik(Nq_CRc);
else
    Nq_CRc.mu=nan;Nq_CRc.sigma=nan;
    NLLq_CRc = nan;
end

%% Find AIC & BIC of fits:
kLin = 2; kQuad = 3;
aicCD1 = 2*NLL_CD + 2*(kLin+1);           % CD
aicCD2 = 2*NLL_CDc + 2*(kLin+1);
aicCD3 = 2*NLLq_CD + 2*(kQuad+1);
aicCD4 = 2*NLLq_CDc + 2*(kQuad+1);
bicCD1 = 2*NLL_CD + log(nSS)*(kLin+1);
bicCD2 = 2*NLL_CDc + log(nSS)*(kLin+1);
bicCD3 = 2*NLLq_CD + log(nSS)*(kQuad+1);
bicCD4 = 2*NLLq_CDc + log(nSS)*(kQuad+1);
% Determine best fit function for CD data:
[~, ix] = min([aicCD1 aicCD2 aicCD3 aicCD4]);
results.CD.bestFitA = ix;
[~, ix] = min([bicCD1 bicCD2 bicCD3 bicCD4]);
results.CD.bestFitB = ix;

aicCR1 = 2*NLL_CR + 2*(kLin+1);           % CR
aicCR2 = 2*NLL_CRc + 2*(kLin+1);
aicCR3 = 2*NLLq_CR + 2*(kQuad+1);
aicCR4 = 2*NLLq_CRc + 2*(kQuad+1);
bicCR1 = 2*NLL_CR + log(nSS)*(kLin+1);
bicCR2 = 2*NLL_CRc + log(nSS)*(kLin+1);
bicCR3 = 2*NLLq_CR + log(nSS)*(kQuad+1);
bicCR4 = 2*NLLq_CRc + log(nSS)*(kQuad+1);
% Determine best fit function for CR data:
[~, ix] = min([aicCR1 aicCR2 aicCR3 aicCR4]);
results.CR.bestFitA = ix;
[~, ix] = min([bicCR1 bicCR2 bicCR3 bicCR4]);
results.CR.bestFitB = ix;

%% Store results in output structure:
results.initials = res.params.initials;
results.tFreq = res.params.tFreq;
if ~any([res.params.dOffset{:}]>0)
    results.bPhSc = false;
else
    results.bPhSc = true;
end
results.CD.data.x = nSetSize;
results.CD.data.mRT = mRT_CD;
results.CD.data.stdRT = stdRT_CD;
results.CD.data.xCloud = x_CD;
results.CD.data.RTCloud = RT_CD;

results.CD.linFit.slope = mLinCD;
results.CD.linFit.int = bLinCD;
results.CD.linFit.NLL = NLL_CD;
results.CD.linFit.resMean = N_CD.mu;
results.CD.linFit.resSig = N_CD.sigma;
results.CD.linFit.AIC = aicCD1;
results.CD.linFit.BIC = bicCD1;

results.CD.linFitCloud.slope = mLinCDc;
results.CD.linFitCloud.int = bLinCDc;
results.CD.linFitCloud.NLL = NLL_CDc;
results.CD.linFitCloud.resMean = N_CDc.mu;
results.CD.linFitCloud.resSig = N_CDc.sigma;
results.CD.linFitCloud.AIC = aicCD2;
results.CD.linFitCloud.BIC = bicCD2;

results.CD.quadFit.x2 = m2QuadCD;
results.CD.quadFit.x = mQuadCD;
results.CD.quadFit.b = bQuadCD;
results.CD.quadFit.NLL = NLLq_CD;
results.CD.quadFit.resMean = Nq_CD.mu;
results.CD.quadFit.resSig = Nq_CD.sigma;
results.CD.quadFit.AIC = aicCD3;
results.CD.quadFit.BIC = bicCD3;

results.CD.quadFitCloud.x2 = m2QuadCDc;
results.CD.quadFitCloud.x = mQuadCDc;
results.CD.quadFitCloud.b = bQuadCDc;
results.CD.quadFitCloud.NLL = NLLq_CDc;
results.CD.quadFitCloud.resMean = Nq_CDc.mu;
results.CD.quadFitCloud.resSig = Nq_CDc.sigma;
results.CD.quadFitCloud.AIC = aicCD4;
results.CD.quadFitCloud.BIC = bicCD4;

results.CR.data.x = nSetSize;
results.CR.data.mRT = mRT_CR;
results.CR.data.stdRT = stdRT_CR;
results.CR.data.xCloud = x_CR;
results.CR.data.RTCloud = RT_CR;

results.CR.linFit.slope = mLinCR;
results.CR.linFit.int = bLinCR;
results.CR.linFit.NLL = NLL_CR;
results.CR.linFit.resMean = N_CR.mu;
results.CR.linFit.resSig = N_CR.sigma;
results.CR.linFit.AIC = aicCR1;
results.CR.linFit.BIC = bicCR1;

results.CR.linFitCloud.slope = mLinCRc;
results.CR.linFitCloud.int = bLinCRc;
results.CR.linFitCloud.NLL = NLL_CRc;
results.CR.linFitCloud.resMean = N_CRc.mu;
results.CR.linFitCloud.resSig = N_CRc.sigma;
results.CR.linFitCloud.AIC = aicCR2;
results.CR.linFitCloud.BIC = bicCR2;

results.CR.quadFit.x2 = m2QuadCR;
results.CR.quadFit.x = mQuadCR;
results.CR.quadFit.b = bQuadCR;
results.CR.quadFit.NLL = NLLq_CR;
results.CR.quadFit.resMean = Nq_CR.mu;
results.CR.quadFit.resSig = Nq_CR.sigma;
results.CR.quadFit.AIC = aicCR3;
results.CR.quadFit.BIC = bicCR3;

results.CR.quadFitCloud.x2 = m2QuadCRc;
results.CR.quadFitCloud.x = mQuadCRc;
results.CR.quadFitCloud.b = bQuadCRc;
results.CR.quadFitCloud.NLL = NLLq_CRc;
results.CR.quadFitCloud.resMean = Nq_CRc.mu;
results.CR.quadFitCloud.resSig = Nq_CRc.sigma;
results.CR.quadFitCloud.AIC = aicCR4;
results.CR.quadFitCloud.BIC = bicCR4;


%% Find correlation coefficient between fit & data:
% Find predicted responses:
yfit = polyval([mLinCD bLinCD],results.CD.data.x);
rtemp = corrcoef(results.CD.data.mRT,yfit);
results.CD.linFit.R = rtemp(2);

yfit = polyval([mLinCDc bLinCDc],results.CD.data.xCloud);
yfull = [];
for cc = nSetSize
    yfull = [yfull results.CD.data.RTCloud{cc}];
end
rtemp = corrcoef(yfull,yfit);
results.CD.linFitCloud.R = rtemp(2);

yfit = polyval([mLinCR bLinCR],results.CR.data.x);
rtemp = corrcoef(results.CR.data.mRT,yfit);
results.CR.linFit.R = rtemp(2);

yfit = polyval([mLinCRc bLinCRc],results.CR.data.xCloud);
yfull = [];
for cc = nSetSize
    yfull = [yfull results.CR.data.RTCloud{cc}];
end
rtemp = corrcoef(yfull,yfit);
results.CR.linFitCloud.R = rtemp(2);

end