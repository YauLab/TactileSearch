function ref = MakeTermReference

ref.C = 1;
ref.linear = 2:17;
ref.TD_withinLimb = [26 39 55 66 80 89 101 108];
ref.DD_withinLimb = [110 123 132 137];
ref.TD_homoBPs = [29 44 58 71 75 86 96 105];
ref.DD_homoBPs = [113 120 126 131];
ref.TD_sameSide = [27 28 41 42 52 53 64 65 81 82 91 92 98 99 106 107];
ref.DD_sameSide = [111 112 117 118 133 134 135 136];
ref.TD_diffSide = [30 31 32 43 45 46 56 57 59 68 69 70 76 77 78 85 87 88 94 95 97 102 103 104];
ref.DD_diffSide = [114 115 116 119 121 122 124 125 127 128 129 130];

% Specify total number of linear & interaction terms:
ref.NLTerms = sort([ref.TD_withinLimb ref.DD_withinLimb ref.TD_homoBPs ...
    ref.DD_homoBPs ref.TD_sameSide ref.DD_sameSide ref.TD_diffSide ref.DD_diffSide]);
ref.nLTerms = length(ref.linear);
ref.nNLTerms = length(ref.NLTerms);

% Also store table indicating which sites are part of each term:
[~,NLSites] = Make2ndPBSDict;
ref.siteIDs = nan([ref.nNLTerms 2]);
for ii = 1:ref.nNLTerms
    ref.siteIDs(ii,:) = NLSites(ref.NLTerms(ii)-17,:)-1;
end

% Split into groups based on how many terms in each:
ref.g = cell([1 8]);
ref.g{1} = ref.TD_diffSide;
ref.g{2} = ref.TD_sameSide;
ref.g{3} = ref.DD_diffSide;
ref.g{4} = ref.DD_sameSide;
ref.g{5} = ref.TD_withinLimb;
ref.g{6} = ref.DD_withinLimb;
ref.g{7} = ref.TD_homoBPs;
ref.g{8} = ref.DD_homoBPs;

ref.gLabels = {'T-D Diff Sides','T-D Same Side','D-D Diff Sides','D-D Same Side',...
    'T-D Within-limb','D-D Within-limb','T-D Homologous BPs','D-D Homologous BPs'};

ref.sortTerms = [ref.g{1} ref.g{2} ref.g{3} ref.g{4} ref.g{5} ref.g{6} ref.g{7} ref.g{8}];

ref.sortTerms2 = [ref.g{2} ref.g{5} ref.g{4} ref.g{6} ref.g{1} ref.g{7} ref.g{3} ref.g{8}];
ref.sortInds2 = [2 5 4 6 1 7 3 8];
ref.gLabels2 = ref.gLabels(ref.sortInds2);

ref.sortTerms3 = [ref.g{2} ref.g{4} ref.g{5} ref.g{6} ref.g{1} ref.g{3} ref.g{7} ref.g{8}];
ref.sortInds3 = [2 4 5 6 1 3 7 8];
ref.gLabels3 = ref.gLabels(ref.sortInds3);

ref.sortSiteIDs = nan([ref.nNLTerms 2]);
for ii = 1:ref.nNLTerms
    ref.sortSiteIDs(ii,:) = NLSites(ref.sortTerms(ii)-17,:)-1;
end


ref.sortSiteIDs2 = nan([ref.nNLTerms 2]);
for ii = 1:ref.nNLTerms
    ref.sortSiteIDs2(ii,:) = NLSites(ref.sortTerms2(ii)-17,:)-1;
end

end