function TwoSiteModelResults(res10,res30,res10ps,res30ps)
% TwoSiteModelResults(res10,res30,res10ps,res30ps)
%
% ejh 04/09/19

%% Plot OMP terms as binary arrays & count histograms:
PlotOMPBinary(res10,res30,res10ps,res30ps);

%% Gather average number of terms across subjects:
nSubj = size(res10,1);
nTerms10 = []; nTerms30 = []; nTerms10ps = []; nTerms30ps = [];
rr10 = []; rr10ps = []; rr30 = []; rr30ps = [];
for ii = 1:nSubj
    t = res10(ii).optBIC.nTerms;
    t = t - 17; % remove 1-site terms
    nTerms10 = [nTerms10 t];
    rr10 = [rr10 res10(ii).optBIC.RMSE];
    
    t = res30(ii).optBIC.nTerms;
    t = t - 17;
    nTerms30 = [nTerms30 t]; %#ok<*AGROW>
    rr30 = [rr30 res30(ii).optBIC.RMSE];

    t = res10ps(ii).optBIC.nTerms;
    t = t - 17;
    nTerms10ps = [nTerms10ps t];
    rr10ps = [rr10ps res10ps(ii).optBIC.RMSE];

    t = res30ps(ii).optBIC.nTerms;
    t = t - 17;
    nTerms30ps = [nTerms30ps t];
    rr30ps = [rr30ps res30ps(ii).optBIC.RMSE];
end

% Mean & sem for number of terms across all sessions:
ntt = [nTerms10(:);nTerms10ps(:);nTerms30(:);nTerms30ps(:)];
disp('Number of terms, across sessions (mean & sem):');
mean(ntt)
std(ntt)/sqrt(length(ntt))

%% Bar plot of # of terms over condition, w/ points for each subj results:
[p1,tbl1,stats1] = anova2([[nTerms10';nTerms10ps'] [nTerms30';nTerms30ps']],nSubj);

figure('color',[1 1 1],'pos',[0 0 600 600]); mSize = 10;
nSubj1 = 12;
gap = 0.05; mh = 0.1; mw = 0.12;
subtightplot(1,1,1,gap,mh,mw); hold on;
bColor = [1 1 1];
bbAll = [nTerms10;nTerms10ps;nTerms30;nTerms30ps];
box off; set(gca,'FontSize',18,'TickDir','out','YGrid','on');
ylabel('Number of Interaction Terms');
bb = nanmean(bbAll,2);
ss = nanstd(bbAll,[],2)/sqrt(nSubj);
bar(bb,'FaceColor',bColor,'BarWidth',0.9);
% Add data points:
xS = -.15:.3/nSubj1:.15; 
xScatter = xS(xS~=0); xScatter = shuffleVector(xScatter,1);
for ii = 1:4
    for jj = 1:nSubj
        plot(xScatter(jj)+ii,bbAll(ii,jj),'ko','MarkerSize',mSize);
    end
end

errorbar(1:4,bb,ss,'k.','LineWidth',4,'CapSize',0);
xticks(1:4);
condStrs = {'10T,sync','10T,phsc','30T,sync','30T,phsc'};
xticklabels(condStrs);
maxY = ceil(nanmax(bbAll(:))/10)*10;
yticks(0:10:maxY); 
axis([0.5 4.5 0 maxY]);

%% % Average over subjects too:
[p2,tbl2,stats2] = anova2([[rr10';rr10ps'] [rr30';rr30ps']],nSubj);
avgRMSE = [rr10' rr10ps' rr30' rr30ps'];
mmR2 = nanmean(avgRMSE,1);
seR2 = nanstd(avgRMSE,[],1)/sqrt(nSubj);
mSize = 10;
figure('color',[1 1 1],'pos',[0 0 600 600]); hold on;
bar(mmR2,'BarWidth',0.9,'FaceColor',[1 1 1]);
% Add subject data points:
xS = -.15:.3/nSubj:.15; 
xScatter = xS(xS~=0); xScatter = shuffleVector(xScatter,1);
for ii = 1:4
    for jj = 1:nSubj
        plot(xScatter(jj)+ii,avgRMSE(jj,ii),'ko','MarkerSize',mSize);
    end
end
errorbar(1:4,mmR2,seR2,'k.','LineWidth',4,'CapSize',0);
box off; set(gca,'FontSize',18,'TickDir','out','YGrid','on'); 
axis([0.5 4.5 0 1.8]);
yticks(0:0.2:1.8);
condStrs = {'T10,sync','T10,phsc','T30,sync','T30,phsc'};
xticks(1:4); xticklabels(condStrs);
ylabel('RMSE');
title('Two-site model');


%% BIC & RMSE over number of terms:
rr = {res10,res10ps,res30,res30ps};
nCond = length(rr); nTerms = length(res10(1).optBIC.coefs);
nSess = nCond * nSubj;
BIClines = nan([nSess nTerms]);
RMSElines = nan([nSess nTerms]);
linRMSE = nan([nSubj nCond]);
ict = 1;
for ii = 1:nSubj
    for jj = 1:nCond
        BIClines(ict,:) = rr{jj}(ii).BIC;
        RMSElines(ict,:) = rr{jj}(ii).RMSE;
        linRMSE(ii,jj) = rr{jj}(ii).linear.RMSE;
        ict = ict + 1;
    end
end

% Average across sessions:
BIC1 = nanmean(BIClines,1); 
BICse = nanstd(BIClines,[],1)/sqrt(nSess);
RMSE1 = nanmean(RMSElines,1);
RMSEse = nanstd(RMSElines,[],1)/sqrt(nSess);

%% Plot
figure('color',[1 1 1],'pos',[0 0 700 600]);
alpha = 0.2; cTA = [0 0 1]; cTP = [1 0 0];
buff = 0.5; xx = -15;
yyaxis right;
stdshade(RMSE1(18:end),RMSEse(18:end),alpha,cTP); hold on;
plot(RMSE1(18:end),'r-'); ylabel('RMSE');

rmse17 = mean(RMSElines(:,17)); rmse17s = std(RMSElines(:,17))/sqrt(nSess);
plot(xx-buff,rmse17,'ro','MarkerFaceColor','r');
errorbar(xx-buff,rmse17,rmse17s,'r.','CapSize',0,'LineWidth',1);

yyaxis left; ylabel('BIC');
stdshade(BIC1(18:end),BICse(18:end),alpha,cTA); hold on;
plot(BIC1(18:end),'b-'); 

bic17 = mean(BIClines(:,17)); bic17s = std(BIClines(:,17))/sqrt(nSess);
plot(xx+buff,bic17,'bo','MarkerFaceColor','b'); 
errorbar(xx+buff,bic17,bic17s,'b.','CapSize',0,'LineWidth',1);

ax = gca; minY = ax.YLim(1); maxY = ax.YLim(2);
yticks(minY:600:maxY);

set(gca,'FontSize',20,'TickDir','out');  box off;

xticks([xx 1 20:20:80]); 
xticklabels({'1-site model','1','20','40','60','80'});
xlim([xx*2 85]);
xlabel('Number of two-site terms')


%% One-site RMSE bar plot:
mmR2 = nanmean(linRMSE,1);
seR2 = nanstd(linRMSE,[],1)/sqrt(nSubj);
[p3,tbl3,stats3] = anova1(linRMSE);

% Mean & sem over all sessions:
mmm = nanmean(linRMSE(:));
sss = nanstd(linRMSE(:))/sqrt(length(linRMSE(:)));

mSize = 10;
figure('color',[1 1 1],'pos',[0 0 600 600]); hold on;
bar(mmR2,'BarWidth',0.9,'FaceColor',[1 1 1]);
% Add subject data points:
xS = -.15:.3/nSubj:.15; 
xScatter = xS(xS~=0); xScatter = shuffleVector(xScatter,1);
for ii = 1:4
    for jj = 1:nSubj
        plot(xScatter(jj)+ii,linRMSE(jj,ii),'ko','MarkerSize',mSize);
    end
end
errorbar(1:4,mmR2,seR2,'k.','LineWidth',4,'CapSize',0);
box off; set(gca,'FontSize',18,'TickDir','out','YGrid','on'); 
axis([0.5 4.5 0 1.8]);
yticks(0:0.2:1.8);
condStrs = {'T10,sync','T10,phsc','T30,sync','T30,phsc'};
xticks(1:4); xticklabels(condStrs);
ylabel('RMSE');
title('One-site model');

end