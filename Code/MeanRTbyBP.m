function MeanRTbyBP(resAll)
% MeanRTbyBP(resAll)

nT = 8;
nR = length(resAll);
nSubj = length(resAll{1});
allSubjRTsites1 = nan([nSubj nR nT]);
allSubjRTsites2 = nan([nSubj nR nT]);
for ii = 1:nSubj
    rtSites1 = nan([nR nT]);
    rtSites2 = nan([nR nT]);

    for jj = 1:nR
        res = resAll{jj}{ii};
        bCorrect = res.response.correct';
        RT1 = res.response.RT;
        % z-score RT:
        RT = (RT1-nanmean(RT1(bCorrect)))/nanstd(RT1(bCorrect));
        bNT = res.params.iTarget == 0;
        SS = res.params.ND; 
        bTP = res.params.iTarget > 0;
        SS(bTP) = SS(bTP)+1;
        b1SS = SS == 1;
        i2 = find(bNT & b1SS & bCorrect);
        ds2 = cell2mat(res.params.dSites(i2));
        for kk = 1:nT
            bT = res.params.iTarget == kk;
            RT1 = RT(bT & bCorrect);
            rtSites1(jj,kk) = nanmean(RT1);
            % Find all 1ND trials with site at this location:
            ix2 = i2(ds2==kk);
            RT2 = RT(ix2);
            rtSites2(jj,kk) = nanmean(RT2);
        end
    end
    allSubjRTsites1(ii,:,:) = rtSites1;
    allSubjRTsites2(ii,:,:) = rtSites2;
end

% Target trials:
mmRT = squeeze(mean(allSubjRTsites1,2));
[pT,tblT,statsT] = anova1(mmRT);


% Distractor trials:
mmRT = squeeze(mean(allSubjRTsites2,2));
[pD,tblD,statsD] = anova1(mmRT);


end