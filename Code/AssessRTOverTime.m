function [p,tbl,stats,allRT] = AssessRTOverTime(res10,res10ps,res30,res30ps,homedir)
% AssessRTOverTime
%
% Find average RT for each session in each subject & run 1-way ANOVA on RT
% over session. (12 subj x 8 sessions) --> disregarding condition/block
%
% ejh 11/01/19


b10first = [1 0 1 0 0 1 0 0 1 1 1 0];
bPhScFirst = [0 0 0 0 0 0 1 1 1 1 1 1];
nSess = 8;
nSubj = length(res10);
allRT = nan([nSubj nSess]);

subjIDs = {'PBS01','PBS02','PBS05','PBS07','PBS08','PBS10',...
    'PBS11','PBS12','PBS13','PBS15','PBS16','PBS17'};
curdir = pwd;

datadir = fullfile(homedir,'Data');

for ii = 1:nSubj
    
    subjdir = fullfile(datadir,subjIDs{ii});
    cd(subjdir);
    subjRTs = nan([1 nSess]);
    
    % Load schedules:
    if b10first(ii)
        if bPhScFirst(ii)
            f1 = 'RunSchedule_Jitter_t10'; r1 = res10ps{ii};
            f2 = 'RunSchedule_Jitter_t30'; r2 = res30ps{ii};
            f3 = 'RunSchedule_t10'; r3 = res10{ii};
            f4 = 'RunSchedule_t30'; r4 = res30{ii};
        else
            f1 = 'RunSchedule_t10'; r1 = res10{ii};
            f2 = 'RunSchedule_t30'; r2 = res30{ii};
            f3 = 'RunSchedule_Jitter_t10'; r3 = res10ps{ii};
            f4 = 'RunSchedule_Jitter_t30'; r4 = res30ps{ii};
        end
    elseif bPhScFirst(ii)
        f1 = 'RunSchedule_Jitter_t30'; r1 = res30ps{ii};
        f2 = 'RunSchedule_Jitter_t10'; r2 = res10ps{ii};
        f3 = 'RunSchedule_t30'; r3 = res30{ii};
        f4 = 'RunSchedule_t10'; r4 = res10{ii};
    else
        f1 = 'RunSchedule_t30'; r1 = res30{ii};
        f2 = 'RunSchedule_t10'; r2 = res10{ii};
        f3 = 'RunSchedule_Jitter_t30'; r3 = res30ps{ii};
        f4 = 'RunSchedule_Jitter_t10'; r4 = res10ps{ii};
    end
    
    % Store correct session order:
    load(f1,'sched');
    ix2 = Find2ndSess(sched);
    subjRTs(1) = FindAvgRT(r1,ix2,0);
    subjRTs(3) = FindAvgRT(r1,ix2,1);

    load(f2,'sched');
    ix2 = Find2ndSess(sched);
    subjRTs(2) = FindAvgRT(r2,ix2,0);
    subjRTs(4) = FindAvgRT(r2,ix2,1); 
    
    load(f3,'sched');
    ix2 = Find2ndSess(sched);
    subjRTs(5) = FindAvgRT(r3,ix2,0);
    subjRTs(7) = FindAvgRT(r3,ix2,1);
    
    load(f4,'sched');
    ix2 = Find2ndSess(sched);
    subjRTs(6) = FindAvgRT(r4,ix2,0);
    subjRTs(8) = FindAvgRT(r4,ix2,1);

    allRT(ii,:) = subjRTs;
    
    cd(curdir);
end

[p,tbl,stats] = anova1(allRT);


end

function RT = FindAvgRT(res,ix,b2nd)

% correct only:
RT1 = res.response.RT;
bCorr = res.response.correct;
RT1(~bCorr) = nan;

% choose which session:
if b2nd
    RT2 = RT1(ix:end);
else
    RT2 = RT1(1:(ix-1));
end

% find average:
RT = nanmean(RT2);

end



function ix2 = Find2ndSess(sched)
% Find first trial index of 2nd session:

% Find trial count of 1st session:
trialCt = 0;
for ii = 1:6
    trialCt = trialCt + length(sched(ii).iTarget);
end
ix2 = trialCt + 1;
end