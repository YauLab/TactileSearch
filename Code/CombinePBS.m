function resAll = CombinePBS(res,bNorm)

if ~exist('bNorm','var'), bNorm = false; end

% Initialize:
resAll = res{1};
nRes = length(res); 
nTotTrials = 0;
for ii = 1:nRes
    nTotTrials = nTotTrials + res{ii}.params.nTrials;
end
resAll.params.nTrials = nTotTrials;
resAll.params.dSites = cell([nTotTrials 1]);
resAll.params.dSites(1:res{1}.params.nTrials) = res{1}.params.dSites;
resAll.params.dOffset = cell([nTotTrials 1]);
resAll.params.dOffset(1:res{1}.params.nTrials) = res{1}.params.dOffset;
resAll.params.sessNum = ones([res{1}.params.nTrials 1]); %
if bNorm
    % Normalize:
    resAll.response.RT = resAll.response.RT/nanmean(resAll.response.RT(resAll.response.correct));
    % Standardize:   (normalize??)
%     resAll.response.RT = (resAll.response.RT-nanmean(resAll.response.RT))/nanstd(resAll.response.RT);
    %resAll.response.RT/nanmean(resAll.response.RT);
end

% Add new field for when sessions are combined across condition:
% resAll.params.tFreq1 = repmat(res{1}.params.tFreq,[1 res{1}.params.nTrials]);

% Add all stim params & responses to output structure:
iCt = res{1}.params.nTrials + 1;
for ii = 2:nRes
    nTrials = res{ii}.params.nTrials;
    resAll.params.iTarget = [resAll.params.iTarget; res{ii}.params.iTarget];
    resAll.params.ND = [resAll.params.ND; res{ii}.params.ND];
    resAll.params.ITI = [resAll.params.ITI; res{ii}.params.ITI];
    resAll.params.dSites(iCt:iCt+nTrials-1) = res{ii}.params.dSites;
    resAll.params.dOffset(iCt:iCt+nTrials-1) = res{ii}.params.dOffset;
    resAll.params.tAmps = [resAll.params.tAmps res{ii}.params.tAmps];
    resAll.params.dAmps = [resAll.params.dAmps res{ii}.params.dAmps];
    resAll.params.sessNum = [resAll.params.sessNum; repmat(ii,[nTrials 1])]; %
    resAll.response.correct = [resAll.response.correct res{ii}.response.correct];
    if bNorm
        RT = res{ii}.response.RT;
        meanRT = nanmean(RT(res{ii}.response.correct));
        resAll.response.RT = [resAll.response.RT RT/meanRT];
    else
        resAll.response.RT = [resAll.response.RT res{ii}.response.RT];
    end
    iCt = iCt+nTrials;
    % Add new field for when sessions are combined across condition:
%     resAll.params.tFreq1 = [resAll.params.tFreq1 repmat(res{ii}.params.tFreq,[1 res{ii}.params.nTrials])];
end
end