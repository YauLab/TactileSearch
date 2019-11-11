function ExportPBS(res10all,res30all,res10ps,res30ps,homedir)
% ExportPBS(res10all,res30all,res10ps,res30ps,homedir)
%
% Collate all data & compile into 1 large dataset. 
% Export out as spreadsheet. (w/ labels)


%% Find total number of correct trials over all subjects:
nTrials = 0; nTrialsJ = 0;
nSubj = length(res10all);
for ii = 1:nSubj
    res = res10all{ii}; if ~isempty(res), nTrials = nTrials + sum(res.response.correct); end
    res = res30all{ii}; if ~isempty(res), nTrials = nTrials + sum(res.response.correct); end
    res = res10ps{ii}; if ~isempty(res), nTrialsJ = nTrialsJ + sum(res.response.correct); end
    res = res30ps{ii}; if ~isempty(res), nTrialsJ = nTrialsJ + sum(res.response.correct); end
end

%% Collect all run schedules:
datadir = fullfile(homedir,'Data');
inits = {'PBS01','PBS02','PBS05','PBS07','PBS08','PBS10','PBS11','PBS12',...
    'PBS13','PBS15','PBS16','PBS17'};
schedName10 = 'RunSchedule_t10.mat';
schedName30 = 'RunSchedule_t30.mat';
schedName10j = 'RunSchedule_Jitter_t10.mat';
schedName30j = 'RunSchedule_Jitter_t30.mat';
orderType = [1 2 1 2 2 1 2 1 2 1 2 1];
tPedalLeft = [1 2 2 1 2 1 2 2 1 1 1 2];
orderNames = {'D30First','D10First'};
pedalNames = {'TPedalLeft','TPedalRight'};
nRuns = 12;
stIdx10 = nan([nSubj nRuns]); stIdx30 = nan([nSubj nRuns]);
stIdx10(:,1) = 1; stIdx30(:,1) = 1;
stIdx10j = nan([nSubj nRuns]); stIdx30j = nan([nSubj nRuns]);
stIdx10j(:,1) = 1; stIdx30j(:,1) = 1;
for ii = 1:nSubj
    
    % Load schedule for T10 runs:
    schedfile = fullfile(datadir,inits{ii},schedName10);
    s10 = load(schedfile,'sched');
    % Load schedule for T30 runs:
    schedfile = fullfile(datadir,inits{ii},schedName30);
    s30 = load(schedfile,'sched');
    
    % Load schedule for T10j runs:
    schedfile = fullfile(datadir,inits{ii},schedName10j);
    s10j = load(schedfile,'sched');
    % Load schedule for T30j runs:
    schedfile = fullfile(datadir,inits{ii},schedName30j);
    s30j = load(schedfile,'sched');
    
    % Find start indices for each run:
    for jj = 1:nRuns-1
        nTrials10 = length(s10.sched(jj).iTarget);
        stIdx10(ii,jj+1) = stIdx10(ii,jj) + nTrials10;
        
        nTrials30 = length(s30.sched(jj).iTarget);
        stIdx30(ii,jj+1) = stIdx30(ii,jj) + nTrials30;
        
        nTrials10j = length(s10j.sched(jj).iTarget);
        stIdx10j(ii,jj+1) = stIdx10j(ii,jj) + nTrials10j;
        
        nTrials30j = length(s30j.sched(jj).iTarget);
        stIdx30j(ii,jj+1) = stIdx30j(ii,jj) + nTrials30j;
        
    end    
end

%% Combine data:
siteNames = {'TA','LE','LI','LK','LF','RE','RI','RK','RF'};
subjIDs = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12'};
iSubj = 1; iRT = 2; iSS = 3; iTP = 4;
iTF = 5; iTS = 6; iTrialNum = 7; iRunNum = 8; iSessNum = 9;
iND = 10; iTA = 11; iSync = 12; iTFO = 13; iTPL = 14; iGrp = 15;

% Add labels:
labels = {'Subject','RT','SetSize','TPresent','TFreq','TSite',...
    'Trial','Run','Session','ND','TAbsent','DSync','TFOrder','TPedal','DSOrder'};
nParams = size(labels,2);
output1 = cell([nTrials+1 nParams]); output1(1,:) = labels;
outputJ = cell([nTrialsJ+1 nParams]); outputJ(1,:) = labels;

trialCt = 2; trialCtJ = 2;
for ii = 1:nSubj
    
    r1 = res10all{ii}; r2 = res10ps{ii};
    r3 = res30all{ii}; r4 = res30ps{ii};
    d1 = str2double(r1.params.date);
    d2 = str2double(r2.params.date);
    d3 = str2double(r3.params.date);
    if d1 > d3
        b10First = false;
    else
        b10First = true;
    end
    if d1 > d2
        bSyncFirst = false;
    else
        bSyncFirst = true;
    end
    
    for jj = 1:2
        if jj == 1
            res = r1;
            resJ = r2;
            tFreq = 'Distractors30Hz';
            sIdx = stIdx10(ii,:);
            sIdxJ = stIdx10j(ii,:);
        else
            res = r3;
            resJ = r4;
            tFreq = 'Distractors10Hz';
            sIdx = stIdx30(ii,:);
            sIdxJ = stIdx30j(ii,:);
        end
        
        if ~isempty(res)
            bCorr = res.response.correct; bCorrJ = resJ.response.correct;
            nKeep = sum(bCorr); nKeepJ = sum(bCorrJ);
            iTarget = res.params.iTarget; iTargetJ = resJ.params.iTarget;
            
            tPresent = iTarget>0;
            tPKeep = tPresent(bCorr); 
            tPStrs = cell(size(tPKeep)); 
            tPStrs(tPKeep) = {'YES'}; tPStrs(~tPKeep) = {'NO'};
            tAbsent = ~tPresent;
            tAKeep = tAbsent(bCorr);
            tAStrs = cell(size(tAKeep)); 
            tAStrs(tAKeep) = {'YES'}; tAStrs(~tAKeep) = {'NO'};
            tPresentJ = iTargetJ>0;
            tPKeepJ = tPresentJ(bCorrJ);
            tPStrsJ = cell(size(tPKeepJ)); 
            tPStrsJ(tPKeepJ) = {'YES'}; tPStrsJ(~tPKeepJ) = {'NO'};
            tAbsentJ = ~tPresentJ;
            tAKeepJ = tAbsentJ(bCorrJ);
            tAStrsJ = cell(size(tAKeepJ)); 
            tAStrsJ(tAKeepJ) = {'YES'}; tAStrsJ(~tAKeepJ) = {'NO'};
            
            RT = res.response.RT(bCorr);
            ND = res.params.ND(bCorr);
            iTKeep = iTarget(bCorr);
            RTJ = resJ.response.RT(bCorrJ);
            NDJ = resJ.params.ND(bCorrJ);
            iTKeepJ = iTargetJ(bCorrJ);
            
            ctInds = 1:length(bCorr);
            ctInds = ctInds(bCorr);
            
            output1(trialCt:trialCt+nKeep-1,iSubj) = subjIDs(ii);
            output1(trialCt:trialCt+nKeep-1,iTF) = {tFreq};
            
            ctIndsJ = 1:length(bCorrJ);
            ctIndsJ = ctIndsJ(bCorrJ);
            
            outputJ(trialCtJ:trialCtJ+nKeepJ-1,iSubj) = subjIDs(ii);
            outputJ(trialCtJ:trialCtJ+nKeepJ-1,iTF) = {tFreq};
            
            for tt = 1:nKeep
                tIdx = trialCt+tt-1;
                output1(tIdx,iRT) = {RT(tt)};
                output1(tIdx,iTP) = tPStrs(tt);
                output1(tIdx,iTA) = tAStrs(tt);
                if tPKeep(tt)
                    output1(tIdx,iSS) = {ND(tt)+1};
                else
                    output1(tIdx,iSS) = {ND(tt)};
                end
                output1(tIdx,iND) = {ND(tt)};
                output1(tIdx,iTS) = siteNames(iTKeep(tt)+1);
                % Determine what run & session this trial is in:
                output1(tIdx,iTrialNum) = {ctInds(tt)};
                iRun = find(ctInds(tt)<sIdx,1,'first') - 1;
                if isempty(iRun), iRun = nRuns; end
                if b10First && bSyncFirst
                    if jj == 1
                        if iRun > 6
                            output1(tIdx,iSessNum) = {3};
                        else
                            output1(tIdx,iSessNum) = {1};
                        end
                    else
                        if iRun > 6
                            output1(tIdx,iSessNum) = {4};
                        else
                            output1(tIdx,iSessNum) = {2};
                        end
                    end
                elseif ~b10First && bSyncFirst
                    if jj == 1
                        if iRun > 6
                            output1(tIdx,iSessNum) = {4};
                        else
                            output1(tIdx,iSessNum) = {2};
                        end
                    else
                        if iRun > 6
                            output1(tIdx,iSessNum) = {3};
                        else
                            output1(tIdx,iSessNum) = {1};
                        end
                    end
                elseif b10First && ~bSyncFirst
                    if jj == 1
                        if iRun > 6
                            output1(tIdx,iSessNum) = {7};
                        else
                            output1(tIdx,iSessNum) = {5};
                        end
                    else
                        if iRun > 6
                            output1(tIdx,iSessNum) = {8};
                        else
                            output1(tIdx,iSessNum) = {6};
                        end
                    end
                elseif ~b10First && ~bSyncFirst
                    if jj == 1
                        if iRun > 6
                            output1(tIdx,iSessNum) = {8};
                        else
                            output1(tIdx,iSessNum) = {6};
                        end
                    else
                        if iRun > 6
                            output1(tIdx,iSessNum) = {7};
                        else
                            output1(tIdx,iSessNum) = {5};
                        end
                    end
                end
                % Make sure all run numbers are w/n-session (1-6):
                iRun1 = iRun; if iRun1 > 6, iRun1 = iRun1 - 6; end
                output1(tIdx,iRunNum) = {iRun1};
                output1(tIdx,iTFO) = orderNames(orderType(ii));
                output1(tIdx,iTPL) = pedalNames(tPedalLeft(ii));
                if ii > 6
                    output1(tIdx,iGrp) = {'PhScFirst'};
                else
                    output1(tIdx,iGrp) = {'SyncFirst'};
                end
            end
            
            for tt = 1:nKeepJ
                tIdx = trialCtJ+tt-1;
                outputJ(tIdx,iRT) = {RTJ(tt)};
                outputJ(tIdx,iTP) = tPStrsJ(tt);
                outputJ(tIdx,iTA) = tAStrsJ(tt);
                if tPKeepJ(tt)
                    outputJ(tIdx,iSS) = {NDJ(tt)+1};
                else
                    outputJ(tIdx,iSS) = {NDJ(tt)};
                end
                outputJ(tIdx,iND) = {NDJ(tt)};
                outputJ(tIdx,iTS) = siteNames(iTKeepJ(tt)+1);
                % Determine what run & session this trial is in:
                outputJ(tIdx,iTrialNum) = {ctIndsJ(tt)};
                iRun = find(ctIndsJ(tt)<sIdxJ,1,'first') - 1;
                if isempty(iRun), iRun = nRuns; end
                if b10First && bSyncFirst
                    if jj == 1
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {7};
                        else
                            outputJ(tIdx,iSessNum) = {5};
                        end
                    else
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {8};
                        else
                            outputJ(tIdx,iSessNum) = {6};
                        end
                    end
                elseif ~b10First && bSyncFirst
                    if jj == 1
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {8};
                        else
                            outputJ(tIdx,iSessNum) = {6};
                        end
                    else
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {7};
                        else
                            outputJ(tIdx,iSessNum) = {5};
                        end
                    end
                elseif b10First && ~bSyncFirst
                    if jj == 1
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {3};
                        else
                            outputJ(tIdx,iSessNum) = {1};
                        end
                    else
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {4};
                        else
                            outputJ(tIdx,iSessNum) = {2};
                        end
                    end
                elseif ~b10First && ~bSyncFirst
                    if jj == 1
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {4};
                        else
                            outputJ(tIdx,iSessNum) = {2};
                        end
                    else
                        if iRun > 6
                            outputJ(tIdx,iSessNum) = {3};
                        else
                            outputJ(tIdx,iSessNum) = {1};
                        end
                    end
                end
                % Make sure all run numbers are w/n-session (1-6):
                iRun1 = iRun; if iRun1 > 6, iRun1 = iRun1 - 6; end
                outputJ(tIdx,iRunNum) = {iRun1};
                outputJ(tIdx,iTFO) = orderNames(orderType(ii));
                outputJ(tIdx,iTPL) = pedalNames(tPedalLeft(ii));
                if ii > 6
                    outputJ(tIdx,iGrp) = {'PhScFirst'};
                else
                    outputJ(tIdx,iGrp) = {'SyncFirst'};
                end
            end
            
            trialCt = trialCt + nKeep;
            trialCtJ = trialCtJ + nKeepJ;
        end
    end
    
end
output1(2:end,iSync) = {'YES'};
outputJ(2:end,iSync) = {'NO'};

%% Export to excel document:
xlswrite('SearchData_Sync.xls',output1);
xlswrite('SearchData_PhSc.xls',outputJ);

end