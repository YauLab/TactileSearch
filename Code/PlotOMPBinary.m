function [meanCoef] = PlotOMPBinary(res10,res30,res10ps,res30ps)

if ~exist('pType','var'), pType = 2; end

nSubj = length(res10);

% Find total number of interaction terms & their indices:
ref = MakeTermReference;
NLTerms = ref.sortTerms2;
grpLabels = ref.gLabels2;
sortInds = ref.sortInds2;

nNLTerms = length(NLTerms);
edges = 1:1:nNLTerms+1;

tmp = false([1 nNLTerms]); tmp1 = nan([1 nNLTerms]);
bins10 = false([nSubj nNLTerms]); bins10ps = false([nSubj nNLTerms]);
bins30 = false([nSubj nNLTerms]); bins30ps = false([nSubj nNLTerms]);
coef10 = nan([nSubj nNLTerms]); coef10ps = nan([nSubj nNLTerms]);
coef30 = nan([nSubj nNLTerms]); coef30ps = nan([nSubj nNLTerms]);
terms10 = []; terms10ps = [];
terms30 = []; terms30ps = [];
for ii = 1:nSubj
    
    t10 = res10(ii).optBIC.terms;
    t30 = res30(ii).optBIC.terms;
    t10ps = res10ps(ii).optBIC.terms;
    t30ps = res30ps(ii).optBIC.terms;
    
    % Make binary vectors indicating presence or absence of each term:
    [tt10,inds] = intersect(NLTerms,t10);
    tmp10 = tmp; c10 = tmp1;
    tmp10(inds) = true;
    c10(inds) = res10(ii).optBIC.coefs(tt10);
    coef10(ii,:) = c10;
    
    [tt30,inds] = intersect(NLTerms,t30);
    tmp30 = tmp; c30 = tmp1;
    tmp30(inds) = true;
    c30(inds) = res30(ii).optBIC.coefs(tt30);
    coef30(ii,:) = c30;
    
    [tt10ps,inds] = intersect(NLTerms,t10ps);
    tmp10ps = tmp; c10ps = tmp1;
    tmp10ps(inds) = true;
    c10ps(inds) = res10ps(ii).optBIC.coefs(tt10ps);
    coef10ps(ii,:) = c10ps;
    
    [tt30ps,inds] = intersect(NLTerms,t30ps);
    tmp30ps = tmp; c30ps = tmp1;
    tmp30ps(inds) = true;
    c30ps(inds) = res30ps(ii).optBIC.coefs(tt30ps);
    coef30ps(ii,:) = c30ps;
    
    % Store in output array for correct condition:
    bins10(ii,:) = tmp10;
    bins10ps(ii,:) = tmp10ps;
    bins30(ii,:) = tmp30;
    bins30ps(ii,:) = tmp30ps;
    
    terms10 = [terms10; t10];
    terms10ps = [terms10ps; t10ps];
    terms30 = [terms30; t30];
    terms30ps = [terms30ps; t30ps]; %#ok<*AGROW>
end

%% Edit arrays so each column (term) is a different color (value):
cmap = cbrewer('qual','Dark2',8); 
cmap = [ 1 1 1; cmap];
toPlot10 = zeros(size(bins10));
toPlot30 = zeros(size(bins30));
toPlot10ps = zeros(size(bins10ps));
toPlot30ps = zeros(size(bins30ps));

% Color code depends on type of interaction term: (8 groups)
nGroups = size(ref.g,2);
ixs = zeros([1 nNLTerms]);
for ii = 1:nNLTerms
    % Find which term this is & which group it's in:
    t = NLTerms(ii);
    for gg = 1:nGroups
        if ~isempty(intersect(t,ref.g{gg}))
            % Find index of this group in group sort inds:
            gi = find(sortInds==gg);
            ix = gi;
        end
    end
    ixs(ii) = ix;
    
    % Color code based on group:
    toPlot10(:,ii) = ix;
    toPlot30(:,ii) = ix;
    toPlot10ps(:,ii) = ix;
    toPlot30ps(:,ii) = ix;
end

%% Turn 10T term indices into counts for histogram:
% Remove all linear terms:
terms10(terms10<=17) = [];
% Turn all terms into their indices in list of NL terms:
inds10 = zeros(size(terms10));
for ii = 1:length(terms10)
    t = terms10(ii);
    inds10(ii) = find(NLTerms == t);
end

% nPltTerms10 = max(inds10);
% [N10,~] = hist(inds10,nPltTerms10);
N10 = histcounts(inds10,edges);

% Remove all linear terms:
terms10ps(terms10ps<=17) = [];
% Turn all terms into their indices in list of NL terms:
inds10ps = zeros(size(terms10ps));
for ii = 1:length(terms10ps)
    t = terms10ps(ii);
    inds10ps(ii) = find(NLTerms == t);
end

% nPltTerms10ps = max(inds10ps);
% [N10ps,~] = hist(inds10ps,nPltTerms10ps);
N10ps = histcounts(inds10ps,edges);

%% Turn 30T term indices into counts for histogram:
% Remove all linear terms:
terms30(terms30<=17) = [];
% Turn all terms into their indices in list of NL terms:
inds30 = zeros(size(terms30));
for ii = 1:length(terms30)
    t = terms30(ii);
    inds30(ii) = find(NLTerms == t);
end

% nPltTerms30 = max(inds30);
% [N30,~] = hist(inds30,nPltTerms30);
N30 = histcounts(inds30,edges);

% Remove all linear terms:
terms30ps(terms30ps<=17) = [];
% Turn all terms into their indices in list of NL terms:
inds30ps = zeros(size(terms30ps));
for ii = 1:length(terms30ps)
    t = terms30ps(ii);
    inds30ps(ii) = find(NLTerms == t);
end

% nPltTerms30ps = max(inds30ps);
% [N30ps,~] = hist(inds30ps,nPltTerms30ps);
N30ps = histcounts(inds30ps,edges);

maxN = max([N10 N30 N10ps N30ps]);

%% Plot w/ subjects next to one another: (interleaved conditions)
cmap2 = cmap(2:end,:);
fSize = 20;

toPlotINTLV = zeros([nSubj*4 nNLTerms]);
toPlotOverAll = zeros([nSubj nNLTerms]);
for ii = 1:nSubj
    y1 = toPlot10(ii,:);
    y2 = toPlot10ps(ii,:);
    y3 = toPlot30(ii,:);
    y4 = toPlot30ps(ii,:);
    
    y = [y1;y2;y3;y4];
    bKeep = any(y,1);
    % Loop through columns to find & store correct values:
    for jj = 1:nNLTerms
        if bKeep(jj)
            x = y(:,jj);
            x = x(x>0);
            toPlotOverAll(ii,jj) = x(1);
        end
    end
    
    toPlotINTLV(ii*4-3,:) = y1;
    toPlotINTLV(ii*4-2,:) = y2;
    toPlotINTLV(ii*4-1,:) = y3;
    toPlotINTLV(ii*4,:) = y4;
    
    coefINTLV(ii*4-3,:) = coef10(ii,:);
    coefINTLV(ii*4-2,:) = coef10ps(ii,:);
    coefINTLV(ii*4-1,:) = coef30(ii,:);
    coefINTLV(ii*4,:) = coef30ps(ii,:);
    
end


%% Find % that each term occurs across sessions:
pOccur = sum(~isnan(coefINTLV))/nNLTerms*100;
disp('Terms occur in what percentage of sessions, m/se:');
mpo = mean(pOccur)
sepo = std(pOccur)/sqrt(nNLTerms) %#ok<*NOPRT,*NASGU>


%% Find mean weights of key terms (split by ipsi vs. contra terms):
% Return vector of top terms (> 1/2 of sessions):
topThresh = size(coefINTLV,1)/2;
indsAll = [inds10;inds10ps;inds30;inds30ps];
Nall = histcounts(indsAll,edges);       % histogram across all sessions
meanCoef = nanmean(coefINTLV,1);

iKeep = find(Nall > topThresh);
ipsi = [ref.TD_withinLimb ref.TD_sameSide ref.DD_withinLimb ref.DD_sameSide];
% contra = [ref.TD_homoBPs ref.TD_diffSide ref.DD_homoBPs ref.DD_diffSide];
ipsiW = []; contraW = [];
for ii = 1:length(iKeep)
    % Find term & determine if ipsi or contra:
    ix = iKeep(ii);
    t = NLTerms(ix);
    if any(ipsi==t)
        % Add weight to ipsi list:
        ipsiW = [ipsiW meanCoef(ix)];
    else
        contraW = [contraW meanCoef(ix)];
    end
end
% Find mean & se of weights:
disp('Mean/se ipsi weights:');
mi = mean(ipsiW)
sei = std(ipsiW)/sqrt(length(ipsiW))
disp('Mean/se contra weights:');
mc = mean(contraW)
sec = std(contraW)/sqrt(length(contraW))
disp('2-sample t-test:')
[~,p,~,stats] = ttest2(ipsiW,contraW);
p
stats

%% Plot interleaved (all subj & sess) in black & white (& gray):

indsAll = [inds10;inds10ps;inds30;inds30ps];
Nall = histcounts(indsAll,edges);       % histogram across all sessions

figure('color',[1 1 1],'pos',[0 0 1200 600]);
gap = .03; mh = 0.1; mw = 0.05;
subtightplot(2,1,2,gap,mh,mw);

cPlot = coefINTLV;
% Make binary:
bb1 = rem(ixs,2)>0;
nSess = size(cPlot,1);
bb1 = repmat(bb1,[nSess 1]);
bb2 = ~bb1; 

bKeep = abs(cPlot)>0;
cPlot(bKeep&bb1) = 1; 
cPlot(bKeep&bb2) = 2;
cPlot(~bKeep) = 0;

imagesc(cPlot);
cw = [1 1 1]; 
c1 = [.55 .55 .55]; c2 = [.3 .3 .3];
colormap([cw;c1;c2]);

hold on;
% Add horizontal lines separating subjects:
xx = 0.5:nNLTerms+0.5;
yy = 0.5:4:nSubj*4+0.75;
for ii = 1:length(yy)
    plot(xx,repmat(yy(ii),size(xx)),'k-');
end
% Add vertical lines separating the term groups:
yy = 0:1:50;
for ii = 2:nGroups
    ti = find(ixs>ii-1,1,'first') - 0.5;
    xx = repmat(ti,size(yy));
    plot(xx,yy,'k-');
end
set(gca,'TickLength',[0 0]);
set(gca,'xtick',[]);
yticks(2.5:4:50); yticklabels({'1','2','3','4','5','6','7','8','9','10','11','12'});
xlabel('Interaction Term'); ylabel('Subjects');
set(gca,'FontSize',fSize);

% Plot interleaved hist (combined indices)
nPltTermsAll = size(toPlotINTLV,2);
X = 1:nPltTermsAll;
subtightplot(2,1,1,gap,mh,mw);
Bh = bar(X,Nall,'FaceColor','flat','BarWidth',1);
for ii = 1:nPltTermsAll
    if rem(ixs(ii),2)>0
        Bh.CData(ii,:) = c1;
    else
        Bh.CData(ii,:) = c2;
    end
end
maxNn = max(Nall(:));
roundY = 2; maxY = ceil(maxNn/roundY)*roundY;
axis([0.5 nNLTerms+0.5 0 maxY]);

box off;
set(gca,'TickDir','out');
set(gca,'TickLength',[0.001 0.001]);
set(gca,'xtick',[]);

yticks(0:10:maxY);
ylabel('Counts'); set(gca,'FontSize',fSize);

end