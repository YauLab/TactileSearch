function res = AvgOverRedundantTrials(res,bDropIncorrect)
% res = AvgOverRedundantTrials(res)
%
% Average RT over all trials that are repeat conditions. (averaging over
% amplitude jitter)
%
% ejh 09/27/18

if ~isempty(res)
    % Drop incorrect trials, if needed:
    if ~exist('bDropIncorrect','var'), bDropIncorrect = true; end
    if bDropIncorrect
        bCorrect = res.response.correct;
        res.params.iTarget = res.params.iTarget(bCorrect);
        res.params.ND = res.params.ND(bCorrect);
        res.params.dSites = res.params.dSites(bCorrect);
        res.params.dOffset = res.params.dOffset(bCorrect);
        res.response.RT = res.response.RT(bCorrect);
        res.response.correct = res.response.correct(bCorrect);
    end
    
    % For each trial, find if there are any other repeats:
    nTrials = length(res.params.iTarget);
    bDone = false([nTrials 1]);
    RT = nan([1 nTrials]);
    nRepeats = [];
    for ii = 1:nTrials
        if ~bDone(ii)
            % Find trial info:
            iT = res.params.iTarget(ii);
            dSites = res.params.dSites{ii};
            
            % Find all trials with this target:
            tIdx = find(res.params.iTarget == iT);
            
            % Loop through trials with this target & find others with these dSites:
            sameInds = [];
            for jj = 1:length(tIdx)
                trial = tIdx(jj);
                if isempty(setxor(dSites,res.params.dSites{trial}))
                    sameInds = [sameInds trial]; %#ok<AGROW>
                end
            end
            
            % Add trial w/ averaged RT to new results structure:
            RT(ii) = nanmean(res.response.RT(sameInds));
            
            nRepeats = [nRepeats length(sameInds)]; %#ok<AGROW>
            
            % Remove trials from further analysis:
            bDone(sameInds) = true;
        end
        
    end
    
    % Remove all trials that have NaN RTs now:
    bNaN = isnan(RT);
    res.params.iTarget = res.params.iTarget(~bNaN);
    res.params.ND = res.params.ND(~bNaN);
    res.params.dSites = res.params.dSites(~bNaN);
    res.params.dOffset = res.params.dOffset(~bNaN);
    res.params.nTrials = length(res.params.iTarget);
    res.response.RT = RT(~bNaN);
    res.response.correct = res.response.correct(~bNaN);
end
end