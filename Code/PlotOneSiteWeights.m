function PlotOneSiteWeights(res10,res30,res10ps,res30ps)
% PlotOneSiteWeights(res10,res30,res10ps,res30ps)


%% Pull one-site term weights (1-17) for each condition:
nSubj = size(res10,1);
c10 = nan([nSubj 17]); c10ps = nan([nSubj 17]);
c30 = nan([nSubj 17]); c30ps = nan([nSubj 17]);
nTerms = 17;
terms = 1:nTerms;
nCond = 4; %    # conditions
nSites = 8; %   # body sites

%% Z-score RTs:
R = nan([nSubj nCond]);
blank = zeros([1 length(res10(1).optBIC.coefs)]);
for ii = 1:nSubj
    
    res = res10(ii).res;
    res.response.RT = (res.response.RT - nanmean(res.response.RT))/nanstd(res.response.RT);
    c10(ii,:) = FindNWeights(res,terms,res10(ii).D);
    c = blank; c(1:17) = c10(ii,:);
    [~, ~, ~, ~, R(ii,1)] = FitNLinReg3(res,c,res10(ii).D,0);
    
    res = res30(ii).res;
    res.response.RT = (res.response.RT - nanmean(res.response.RT))/nanstd(res.response.RT);
    c30(ii,:) = FindNWeights(res,terms,res30(ii).D);
    c = blank; c(1:17) = c30(ii,:);
    [~, ~, ~, ~, R(ii,2)] = FitNLinReg3(res,c,res30(ii).D,0);
    
    res = res10ps(ii).res;
    res.response.RT = (res.response.RT - nanmean(res.response.RT))/nanstd(res.response.RT);
    c10ps(ii,:) = FindNWeights(res,terms,res10ps(ii).D);
    c = blank; c(1:17) = c10ps(ii,:);
    [~, ~, ~, ~, R(ii,3)] = FitNLinReg3(res,c,res10ps(ii).D,0);
    
    res = res30ps(ii).res;
    res.response.RT = (res.response.RT - nanmean(res.response.RT))/nanstd(res.response.RT);
    c30ps(ii,:) = FindNWeights(res,terms,res30ps(ii).D);
    c = blank; c(1:17) = c30ps(ii,:);
    [~, ~, ~, ~, R(ii,4)] = FitNLinReg3(res,c,res30ps(ii).D,0);
end


%% Reorder coefficients, so x-axis (8 terms) of body map is symmetric:
% Target terms:
tTerms = nan([nSubj nSites nCond]);
tTerms(:,:,1) = c10(:,[2:5 9 8 7 6]);
tTerms(:,:,2) = c10ps(:,[2:5 9 8 7 6]);
tTerms(:,:,3) = c30(:,[2:5 9 8 7 6]);
tTerms(:,:,4) = c30ps(:,[2:5 9 8 7 6]);

% Distractor terms:
dTerms = nan([nSubj nSites nCond]);
dTerms(:,:,1) = c10(:,[10:13 17 16 15 14]);
dTerms(:,:,2) = c10ps(:,[10:13 17 16 15 14]);
dTerms(:,:,3) = c30(:,[10:13 17 16 15 14]);
dTerms(:,:,4) = c30ps(:,[10:13 17 16 15 14]);


%% Find mean weights (avgd w/n subject across conditions first):
avgT = nanmean(tTerms,3);
avgD = nanmean(dTerms,3);
modAvgD = avgD(:,[1 3:6 8]);

[pT,tblT,sT] = anova1(avgT);
[pD,tblD,sD] = anova1(avgD);
[pDmod,tblDmod,sDmod] = anova1(modAvgD); %#ok<*ASGLU>

%% Split into separate T & D plots -> show data points!!:
tColor = [0.9 0.5 0.5]; dColor = [0.6 0.6 0.9]; mColor = [0 0 0];
gap = 0.1; mh = 0.1; mw = 0.125;
fSize = 20;

% T plot:
figure('color',[1 1 1],'pos',[0 0 1200 800]);
subtightplot(1,1,1,gap,mh,mw); hold on;
meanTP1 = nanmean(avgT);
stdTP1 = nanstd(avgT);
nSess = size(avgT,1);
seTP1 = stdTP1/sqrt(nSess);
sTP = size(avgT,2);
box off; set(gca,'FontSize',fSize,'TickDir','out','YGrid','on');
ylabel('Weight');
bar(meanTP1,'FaceColor',tColor,'BarWidth',0.9);

% Add data points:
mSize = 10; lWidth = 4;
% xS = -.15:.3/nSess:.15;
xS = -.1:.2/nSubj:.1;
xScatter = xS(xS~=0); xScatter = shuffleVector(xScatter,1);

for ii = 1:8
    for jj = 1:nSubj
        mD = nanmean(tTerms(jj,ii,:),3);
        plot(xScatter(jj)+ii,mD,'ko','MarkerEdgeColor',mColor,'MarkerSize',mSize);
    end
end

errorbar(1:sTP,meanTP1,seTP1,'k.','LineWidth',lWidth,'CapSize',0);
roundY = 0.5;
minY = floor(nanmin(avgT(:))/roundY)*roundY;
maxY = ceil(nanmax(avgT(:))/roundY)*roundY;
axis([0.25 8.75 minY maxY]);
xticks(1:sTP); yticks(-2:1:1);
xticklabels({'LE','LI','LK','LF','RF','RK','RI','RE'});
title('Target Terms','FontSize',fSize);
xlabel('Term');

%% D plot:
figure('color',[1 1 1],'pos',[0 0 1200 800]);
subtightplot(1,1,1,gap,mh,mw); hold on;
meanTP2 = nanmean(avgD);
stdTP2 = nanstd(avgD);
seTP2 = stdTP2/sqrt(nSess);
sTP = size(avgD,2);
box off; set(gca,'FontSize',fSize,'TickDir','out','YGrid','on');
ylabel('Weight');
bar(meanTP2,'FaceColor',dColor,'BarWidth',0.9);

% Add data points:
xS = -.1:.2/nSubj:.1;
xScatter = xS(xS~=0); xScatter = shuffleVector(xScatter,1);
for ii = 1:8
    for jj = 1:nSubj
        mD = nanmean(dTerms(jj,ii,:),3);
        plot(xScatter(jj)+ii,mD,'ko','MarkerEdgeColor',mColor,'MarkerSize',mSize);
    end
end

errorbar(1:sTP,meanTP2,seTP2,'k.','LineWidth',lWidth,'CapSize',0);
roundY = 0.1;
minY = floor(nanmin(avgD(:))/roundY)*roundY;
maxY = ceil(nanmax(avgD(:))/roundY)*roundY;
axis([0.25 8.75 minY maxY]);
xticks(1:sTP); yticks(minY:0.1:maxY);
xticklabels({'LE','LI','LK','LF','RF','RK','RI','RE'});
title('Distractor Terms','FontSize',fSize);
xlabel('Term');

end