function [res10all, res30all] = CollectPhScResults(bAvg,bModify,bDrop,homedir)
if ~exist('bAvg','var'), bAvg = true; end
if ~exist('bModify','var'), bModify = true; end
if ~exist('bDrop','var'), bDrop = true; end

r10_1=CollateSearchData('pbs01',10,30,1,homedir);
r10_2=CollateSearchData('pbs02',10,30,1,homedir);
r10_5=CollateSearchData('pbs05',10,30,1,homedir);
r10_7=CollateSearchData('pbs07',10,30,1,homedir);
r10_8=CollateSearchData('pbs08',10,30,1,homedir);
r10_10=CollateSearchData('pbs10',10,30,1,homedir);
r10_11=CollateSearchData('pbs11',10,30,1,homedir);
r10_12=CollateSearchData('pbs12',10,30,1,homedir);
r10_13=CollateSearchData('pbs13',10,30,1,homedir);
r10_15=CollateSearchData('pbs15',10,30,1,homedir);
r10_16=CollateSearchData('pbs16',10,30,1,homedir);
r10_17=CollateSearchData('pbs17',10,30,1,homedir);

r30_1=CollateSearchData('pbs01',30,10,1,homedir);
r30_2=CollateSearchData('pbs02',30,10,1,homedir);
r30_5=CollateSearchData('pbs05',30,10,1,homedir);
r30_7=CollateSearchData('pbs07',30,10,1,homedir);
r30_8=CollateSearchData('pbs08',30,10,1,homedir);
r30_10=CollateSearchData('pbs10',30,10,1,homedir);
r30_11=CollateSearchData('pbs11',30,10,1,homedir);
r30_12=CollateSearchData('pbs12',30,10,1,homedir);
r30_13=CollateSearchData('pbs13',30,10,1,homedir);
r30_15=CollateSearchData('pbs15',30,10,1,homedir);
r30_16=CollateSearchData('pbs16',30,10,1,homedir);
r30_17=CollateSearchData('pbs17',30,10,1,homedir);

if bAvg
    r10_1=AvgOverRedundantTrials(r10_1,bDrop);
    r10_2=AvgOverRedundantTrials(r10_2,bDrop);
    r10_5=AvgOverRedundantTrials(r10_5,bDrop);
    r10_7=AvgOverRedundantTrials(r10_7,bDrop);
    r10_8=AvgOverRedundantTrials(r10_8,bDrop);
    r10_10=AvgOverRedundantTrials(r10_10,bDrop);
    r10_11=AvgOverRedundantTrials(r10_11,bDrop);
    r10_12=AvgOverRedundantTrials(r10_12,bDrop);
    r10_13=AvgOverRedundantTrials(r10_13,bDrop);
    r10_15=AvgOverRedundantTrials(r10_15,bDrop);
    r10_16=AvgOverRedundantTrials(r10_16,bDrop);
    r10_17=AvgOverRedundantTrials(r10_17,bDrop);

    r30_1=AvgOverRedundantTrials(r30_1,bDrop);
    r30_2=AvgOverRedundantTrials(r30_2,bDrop);
    r30_5=AvgOverRedundantTrials(r30_5,bDrop);
    r30_7=AvgOverRedundantTrials(r30_7,bDrop);
    r30_8=AvgOverRedundantTrials(r30_8,bDrop);
    r30_10=AvgOverRedundantTrials(r30_10,bDrop);
    r30_11=AvgOverRedundantTrials(r30_11,bDrop);
    r30_12=AvgOverRedundantTrials(r30_12,bDrop);
    r30_13=AvgOverRedundantTrials(r30_13,bDrop);
    r30_15=AvgOverRedundantTrials(r30_15,bDrop);
    r30_16=AvgOverRedundantTrials(r30_16,bDrop);
    r30_17=AvgOverRedundantTrials(r30_17,bDrop);
end

if bModify
    r10_1 = ModifyTargetIDs(r10_1);
    r10_2 = ModifyTargetIDs(r10_2);
    r10_5 = ModifyTargetIDs(r10_5);
    r10_7 = ModifyTargetIDs(r10_7);
    r10_8 = ModifyTargetIDs(r10_8);
    r10_10 = ModifyTargetIDs(r10_10);
    r10_11 = ModifyTargetIDs(r10_11);
    r10_12 = ModifyTargetIDs(r10_12);
    r10_13 = ModifyTargetIDs(r10_13);
    r10_15 = ModifyTargetIDs(r10_15);
    r10_16 = ModifyTargetIDs(r10_16);
    r10_17 = ModifyTargetIDs(r10_17);
    
    r30_1 = ModifyTargetIDs(r30_1);
    r30_2 = ModifyTargetIDs(r30_2);  
    r30_5 = ModifyTargetIDs(r30_5);
    r30_7 = ModifyTargetIDs(r30_7);
    r30_8 = ModifyTargetIDs(r30_8);
    r30_10 = ModifyTargetIDs(r30_10);
    r30_11 = ModifyTargetIDs(r30_11);
    r30_12 = ModifyTargetIDs(r30_12);
    r30_13 = ModifyTargetIDs(r30_13);
    r30_15 = ModifyTargetIDs(r30_15);
    r30_16 = ModifyTargetIDs(r30_16);
    r30_17 = ModifyTargetIDs(r30_17);
end

res10all = {r10_1,r10_2,r10_5,r10_7,r10_8,r10_10,r10_11,r10_12,r10_13,r10_15,r10_16,r10_17};
res30all = {r30_1,r30_2,r30_5,r30_7,r30_8,r30_10,r30_11,r30_12,r30_13,r30_15,r30_16,r30_17};

end

