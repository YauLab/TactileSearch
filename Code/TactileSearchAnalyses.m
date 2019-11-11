%% Final top-level analysis for Tactile Search project
%  ejh 11/07/19

%%%%%%%%%%%%% Change this path to run all subsequent code %%%%%%%%%%%%%%%%%
homedir = '/Users/liz/Desktop/TactileSearch';


%% Collect all data:
nSubj = 12;

% Unaveraged data:
bAvg = 0; bModify = 1; bDrop = 0;
[res10all,res30all] = CollectSyncResults(bAvg,bModify,bDrop,homedir);
[res10ps,res30ps] = CollectPhScResults(bAvg,bModify,bDrop,homedir);

% Averaged (for unique patterns) correct-only data:
bAvg2 = true; bModify2 = true; bDrop2 = true;
[res10all2,res30all2] = CollectSyncResults(bAvg2,bModify2,bDrop2,homedir);
[res10ps2,res30ps2] = CollectPhScResults(bAvg2,bModify2,bDrop2,homedir);

%% Find percentage of trials with 10 s responses (timed out):
rr = {res10all,res30all,res10ps,res30ps};
RTall = [];
for ii = 1:4
    r = rr{ii};
    for jj = 1:12
        res = r{jj};
        RTall = [RTall res.response.RT];
    end
end
nTotalTrials = length(RTall);
nTimedOut = sum(isnan(RTall));
pTimedOut = nTimedOut/nTotalTrials*100 %#ok<*NOPTS>


%% Performance:
PlotAllAccuracy(res10all,res30all,res10ps,res30ps); % Figure 2


%% Export data for use in R:
% Linear mixed effects modeling and Figures 3 & 4A
ExportPBS(res10all,res30all,res10ps,res30ps,homedir);


%% 1-sample t-test for change in slope due to scrambling: (Figure 4B)
bToss = true;  % correct trials only
[p4,stats4] = LinFits_deltaM(res10all,res30all,res10ps,res30ps,bToss); 


%% Load one-site & two-site model results:
% These results include the one-site model and the two-site models from 
% the matching pursuit model selection procedure. The OMP procedure takes 
% many hours, so I am loading the results here. The associated OMP code 
% is commented out at the bottom of this cell ("RUN OMP").
r10_1 = load(fullfile(homedir,'sAll_10T_Sync_FixLinOMP_10000iter')); 
r10 = r10_1.results;
r30_1 = load(fullfile(homedir,'sAll_30T_Sync_FixLinOMP_10000iter')); 
r30 = r30_1.results;
r10ps1 = load(fullfile(homedir,'sAll_10T_PhSc_FixLinOMP_10000iter')); 
r10ps = r10ps1.results;
r30ps1 = load(fullfile(homedir,'sAll_30T_PhSc_FixLinOMP_10000iter')); 
r30ps = r30ps1.results;

% Part way through analysis, we decided to employ a new method of choosing 
% the best model, finding BIC separately for target present and target 
% absent trials and minimizing their sum. This code modifies past results 
% to be compliant. Run2siteOMP.m (below) is up-to-date.
r10 = MinBICSplitTPTA(r10);
r30 = MinBICSplitTPTA(r30);
r10ps = MinBICSplitTPTA(r10ps);
r30ps = MinBICSplitTPTA(r30ps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN OMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bFixLin = true; subjIDs = 1:12;
% % Run OMP on 10T, Sync sessions:
% results10sync = Run2siteOMP(res10all2,subjIDs,bFixLin);
% % Run OMP on 30T, Sync sessions:
% results30sync = Run2siteOMP(res30all2,subjIDs,bFixLin);
% % Run OMP on 10T, PhSc sessions:
% results10ps = Run2siteOMP(res10ps2,subjIDs,bFixLin);
% % Run OMP on 30T, PhSc sessions:
% results30ps = Run2siteOMP(res30ps2,subjIDs,bFixLin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% One-site model of RT: 
% 1) 1-way ANOVAs on T & D weights
% 2) Plot target and distractor weights (Figure 5)
PlotOneSiteWeights(r10,r30,r10ps,r30ps);


%% RTs by target site (all set sizes) or distractor site (set size = 1):
% Includes 1-way ANOVAs of RTs by target location and of RTs for
% 1-distractor trials by distractor location
MeanRTbyBP({res10all,res30all,res10ps,res30ps});


%% Plot two-site model results:
% Includes 2-way ANOVA of # of two-site terms & of two-site model RMSEs
% Figures 6 & 7
% Includes specifics about two-site terms & their weights
% Includes 1-way ANOVA of one-site model RMSEs
TwoSiteModelResults(r10,r30,r10ps,r30ps);


%% Hamming distance analyses:
TopHammingOMP(r10,r30,r10ps,r30ps);


%% 1-way ANOVA of mean RT over sessions:
[p,tbl,stats,allRT] = AssessRTOverTime(res10all,res10ps,res30all,res30ps,homedir);

