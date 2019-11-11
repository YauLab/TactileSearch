function [p4,stats4] = LinFits_deltaM(res10all,res30all,res10ps,res30ps,bToss)

if ~exist('bToss','var') || isempty(bToss), bToss = true; end

%% Find best fit linear functions for all datasets:
nSubj = length(res10all);
resSync10 = cell([1 nSubj]); resSync30 = cell([1 nSubj]);
resPhSc10 = cell([1 nSubj]); resPhSc30 = cell([1 nSubj]);
for ii = 1:nSubj
    % 10 Hz Target Sync Results:
    r10 = res10all{ii};
    resSync10{ii} = FitRTvSS(r10,bToss);
    
    % 30 Hz Target Sync Results:
    r30 = res30all{ii};
    resSync30{ii} = FitRTvSS(r30,bToss);
    
    % 10 Hz Target PhSc Results:
    r10ps = res10ps{ii};
    resPhSc10{ii} = FitRTvSS(r10ps,bToss);
    
    % 10 Hz PhSc Results:
    r30ps = res30ps{ii};
    resPhSc30{ii} = FitRTvSS(r30ps,bToss);
end

%% Get linear fits:
nSubj = length(res10all);
m1=nan([1 nSubj]);m2=nan([1 nSubj]);m3=nan([1 nSubj]);m4=nan([1 nSubj]);
m5=nan([1 nSubj]);m6=nan([1 nSubj]);m7=nan([1 nSubj]);m8=nan([1 nSubj]);
b1=nan([1 nSubj]);b2=nan([1 nSubj]);b3=nan([1 nSubj]);b4=nan([1 nSubj]);
b5=nan([1 nSubj]);b6=nan([1 nSubj]);b7=nan([1 nSubj]);b8=nan([1 nSubj]);
mRT1=nan([1 nSubj]);mRT2=nan([1 nSubj]);mRT3=nan([1 nSubj]);mRT4=nan([1 nSubj]);
mRT5=nan([1 nSubj]);mRT6=nan([1 nSubj]);mRT7=nan([1 nSubj]);mRT8=nan([1 nSubj]);
mRT10 = cell([1 nSubj]); mRT30 = cell([1 nSubj]);
rCD10 = nan([1 nSubj]); rCD30 = nan([1 nSubj]);
rCR10 = nan([1 nSubj]); rCR30 = nan([1 nSubj]);
rCD10ps = nan([1 nSubj]); rCD30ps = nan([1 nSubj]);
rCR10ps = nan([1 nSubj]); rCR30ps = nan([1 nSubj]);
for ii = 1:nSubj
    if ~isempty(res10all{ii})
        % Slopes: (10T sync, 10T phsc, 30T sync, 30T phsc)
        if ~isempty(resSync10{ii})
            m1(ii) = resSync10{ii}.CD.linFitCloud.slope;
            m5(ii) = resSync10{ii}.CR.linFitCloud.slope;
            b1(ii) = resSync10{ii}.CD.linFitCloud.int;
            b5(ii) = resSync10{ii}.CR.linFitCloud.int;
            % Also grab mean response time for these condition:
            [mRT1(ii),RT1] = FindMeanRTofCloud(resSync10{ii}.CD.data.RTCloud);
            [mRT5(ii),RT5] = FindMeanRTofCloud(resSync10{ii}.CR.data.RTCloud);
            mRT10{ii} = [mRT10{ii} RT1 RT5];
            % GoF of line:
            rCD10(ii) = resSync10{ii}.CD.linFit.R;
            rCR10(ii) = resSync10{ii}.CR.linFit.R;
        end
        if ~isempty(resPhSc10{ii})
            m2(ii) = resPhSc10{ii}.CD.linFitCloud.slope;
            m6(ii) = resPhSc10{ii}.CR.linFitCloud.slope;
            b2(ii) = resPhSc10{ii}.CD.linFitCloud.int;
            b6(ii) = resPhSc10{ii}.CR.linFitCloud.int;
            % Also grab mean response time for these condition:
            mRT2(ii) = FindMeanRTofCloud(resPhSc10{ii}.CD.data.RTCloud);
            mRT6(ii) = FindMeanRTofCloud(resPhSc10{ii}.CR.data.RTCloud);
            % GoF of line:
            rCD10ps(ii) = resPhSc10{ii}.CD.linFit.R;
            rCR10ps(ii) = resPhSc10{ii}.CR.linFit.R;
        end
        if ~isempty(resSync30{ii})
            m3(ii) = resSync30{ii}.CD.linFitCloud.slope;
            m7(ii) = resSync30{ii}.CR.linFitCloud.slope;
            b3(ii) = resSync30{ii}.CD.linFitCloud.int;
            b7(ii) = resSync30{ii}.CR.linFitCloud.int;
            % Also grab mean response time for these condition:
            [mRT3(ii),RT3] = FindMeanRTofCloud(resSync30{ii}.CD.data.RTCloud);
            [mRT7(ii),RT7] = FindMeanRTofCloud(resSync30{ii}.CR.data.RTCloud);
            mRT30{ii} = [mRT10{ii} RT3 RT7];
            % GoF of line:
            rCD30(ii) = resSync30{ii}.CD.linFit.R;
            rCR30(ii) = resSync30{ii}.CR.linFit.R;
        end
        if ~isempty(resPhSc30{ii})
            m4(ii) = resPhSc30{ii}.CD.linFitCloud.slope;
            m8(ii) = resPhSc30{ii}.CR.linFitCloud.slope;
            b4(ii) = resPhSc30{ii}.CD.linFitCloud.int;
            b8(ii) = resPhSc30{ii}.CR.linFitCloud.int;
            % Also grab mean response time for these condition:
            mRT4(ii) = FindMeanRTofCloud(resPhSc30{ii}.CD.data.RTCloud);
            mRT8(ii) = FindMeanRTofCloud(resPhSc30{ii}.CR.data.RTCloud);
            % GoF of line:
            rCD30ps(ii) = resPhSc30{ii}.CD.linFit.R;
            rCR30ps(ii) = resPhSc30{ii}.CR.linFit.R;
        end
    end
end


%% Bar plot of delta slope in each condition w/ dots for each subject's m:
figure('color',[1 1 1],'pos',[0 0 600 600]); mSize = 7.5;
gap = 0.05; mh = 0.1; mw = 0.1;
subtightplot(1,1,1,gap,mh,mw); hold on;
bb = [m2-m1;m6-m5;m4-m3;m8-m7];
meanMs = nanmean(bb,2);
stdMs = nanstd(bb,1,2);
seMs = stdMs/sqrt(nSubj);
box off; set(gca,'FontSize',18,'TickDir','out','YGrid','on');
ylabel('\Delta Slope');

bColor = [.7 .7 .7];
bar(meanMs,'FaceColor',bColor,'BarWidth',0.9);%,'FaceAlpha',0);

% Add data points:
xS = -.15:.3/nSubj:.15;
xScatter = xS(xS~=0); xScatter = shuffleVector(xScatter,1);
for ii = 1:4
    for jj = 1:nSubj
        cc2 = 'k';%cmap(jj,:);
        plot(xScatter(jj)+ii,bb(ii,jj),'ko','MarkerSize',mSize,...
            'MarkerEdgeColor',cc2);
    end
end
errorbar(1:4,meanMs,seMs,'k.','LineWidth',4,'CapSize',0);
roundY = 0.02;
minY = floor(min(bb(:))/roundY)*roundY;
maxY = ceil(max(bb(:))/roundY)*roundY;
axis([0.5 5-0.5 minY maxY]);
xticks(1:4);
xticklabels({'10,TP','10,TA','30,TP','30,TA'});
supertitle({'Average \Delta slope across subjects';'(Scrambled-Aligned)'},0.92,'FontSize',18);

% Perform 1-sample t-test for each distribution being different than zero:
[~,p1,~,~] = ttest(bb(1,:));
[~,p2,~,~] = ttest(bb(2,:));
[~,p3,~,~] = ttest(bb(3,:));
[~,p4,~,stats4] = ttest(bb(4,:));   % stats has t-stat & deg of freedom 

fprintf('\n 1-sample t-tests on delta slopes in each condition:');
fprintf(['\n\n t(10,TP) \n p = ' num2str(p1)]);
fprintf(['\n\n t(10,TA) \n p = ' num2str(p2)]);
fprintf(['\n\n t(30,TP) \n p = ' num2str(p3)]);
fprintf(['\n\n t(30,TA) \n p = ' num2str(p4) '\n']);

end

function [mRT,RT] = FindMeanRTofCloud(allRT)
RT = [];
for ii = 1:8
    RT = [RT allRT{ii}]; %#ok<AGROW>
end
RT = RT*1000;
mRT = nanmean(RT);
end
