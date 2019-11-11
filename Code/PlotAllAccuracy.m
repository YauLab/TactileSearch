function PlotAllAccuracy(res10,res30,res10ps,res30ps)

if ~iscell(res10)
    res10 = {res10}; res30 = {res30};
    res10ps = {res10ps}; res30ps = {res30ps};
end
nSess = length(res10);
% Pre-allocate:
pCD_10 = nan([1 nSess]); pCR_10 = nan([1 nSess]);
pCD_30 = nan([1 nSess]); pCR_30 = nan([1 nSess]);
pCD_10ps = nan([1 nSess]); pCR_10ps = nan([1 nSess]);
pCD_30ps = nan([1 nSess]); pCR_30ps = nan([1 nSess]);
for ii = 1:nSess
    r10 = res10{ii}; r30 = res30{ii};
    r10ps = res10ps{ii}; r30ps = res30ps{ii};
    bCorr10 = r10.response.correct; bCorr30 = r30.response.correct;
    bCorr10ps = r10ps.response.correct; bCorr30ps = r30ps.response.correct;
    bTarg10 = r10.params.iTarget'>0; bTarg30 = r30.params.iTarget'>0;
    bTarg10ps = r10ps.params.iTarget'>0; bTarg30ps = r30ps.params.iTarget'>0;
    % Find percent correct for correct detect, 10 Hz target:
    pCD_10(ii) = sum(bCorr10&bTarg10)/sum(bTarg10);
    pCR_10(ii) = sum(bCorr10&~bTarg10)/sum(~bTarg10);
    pCD_30(ii) = sum(bCorr30&bTarg30)/sum(bTarg30);
    pCR_30(ii) = sum(bCorr30&~bTarg30)/sum(~bTarg30);
    pCD_10ps(ii) = sum(bCorr10ps&bTarg10ps)/sum(bTarg10ps);
    pCR_10ps(ii) = sum(bCorr10ps&~bTarg10ps)/sum(~bTarg10ps);
    pCD_30ps(ii) = sum(bCorr30ps&bTarg30ps)/sum(bTarg30ps);
    pCR_30ps(ii) = sum(bCorr30ps&~bTarg30ps)/sum(~bTarg30ps);
end

pCD_10 = pCD_10*100;
pCR_10 = pCR_10*100;
pCD_30 = pCD_30*100;
pCR_30 = pCR_30*100;
pCD_10ps = pCD_10ps*100;
pCR_10ps = pCR_10ps*100;
pCD_30ps = pCD_30ps*100;
pCR_30ps = pCR_30ps*100;

%% Plot:
nSubj = 12;
mCD10 = mean(pCD_10); sCD10 = std(pCD_10)/sqrt(nSubj);
mCD30 = mean(pCD_30); sCD30 = std(pCD_30)/sqrt(nSubj);
mCR10 = mean(pCR_10); sCR10 = std(pCR_10)/sqrt(nSubj);
mCR30 = mean(pCR_30); sCR30 = std(pCR_30)/sqrt(nSubj);

mCD10ps = mean(pCD_10ps); sCD10ps = std(pCD_10ps)/sqrt(nSubj);
mCD30ps = mean(pCD_30ps); sCD30ps = std(pCD_30ps)/sqrt(nSubj);
mCR10ps = mean(pCR_10ps); sCR10ps = std(pCR_10ps)/sqrt(nSubj);
mCR30ps = mean(pCR_30ps); sCR30ps = std(pCR_30ps)/sqrt(nSubj);

figure('color',[1 1 1],'pos',[0 0 1200 600]); mSize = 10;
gap = 0.1; mh = 0.1; mw = 0.12;

subtightplot(1,2,1,gap,mh,mw); hold on;
bColor = [1 1 1];
bb = [mCD10 mCR10 mCD30 mCR30];
box off; set(gca,'FontSize',18,'TickDir','out','YGrid','on');
ylabel('Accuracy');
bar(bb,'FaceColor',bColor,'BarWidth',0.9);
% Add data points:
xS = -.15:.30/nSubj:.15;
xScatter = xS(xS~=0); xScatter = shuffleVector(xScatter,1);
% buff = 0.15;
for jj = 1:nSubj
    plot(xScatter(jj)+1,pCD_10(jj),'ko','MarkerSize',mSize);
    plot(xScatter(jj)+2,pCR_10(jj),'ko','MarkerSize',mSize);
    plot(xScatter(jj)+3,pCD_30(jj),'ko','MarkerSize',mSize);
    plot(xScatter(jj)+4,pCR_30(jj),'ko','MarkerSize',mSize);
end
errorbar(1:4,[mCD10 mCR10 mCD30 mCR30],[sCD10 sCR10 sCD30 sCR30],'k.','LineWidth',4,'CapSize',0);
xticks([1.5 3.5]);
condStrs = {'10T','30T'};
xticklabels(condStrs);
yticks(0:20:100); 
axis([0.5 4.5 0 105]); set(gca,'FontSize',24);
title('Phase-aligned');

subtightplot(1,2,2,gap,mh,mw); hold on;
bb = [mCD10ps mCR10ps mCD30ps mCR30ps];
box off; set(gca,'FontSize',18,'TickDir','out','YGrid','on');
ylabel('Accuracy');
bar(bb,'FaceColor',bColor,'BarWidth',0.9);
% Add data points:
for jj = 1:nSubj
    plot(xScatter(jj)+1,pCD_10ps(jj),'ko','MarkerSize',mSize);
    plot(xScatter(jj)+2,pCR_10ps(jj),'ko','MarkerSize',mSize);
    plot(xScatter(jj)+3,pCD_30ps(jj),'ko','MarkerSize',mSize);
    plot(xScatter(jj)+4,pCR_30ps(jj),'ko','MarkerSize',mSize);
end
errorbar(1:4,[mCD10ps mCR10ps mCD30ps mCR30ps],[sCD10ps sCR10ps sCD30ps sCR30ps],'k.','LineWidth',4,'CapSize',0);
xticks([1.5 3.5]);
condStrs = {'10T','30T'};
xticklabels(condStrs);
yticks(0:20:100);
axis([0.5 4.5 0 105]); set(gca,'FontSize',24);
title('Phase-scrambled');

end