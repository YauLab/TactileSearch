function [RMSE, R2, simRTs, f1, R] = FitNLinReg3(res,weights,dict,bPlot,nTerms,...
    nLinearTerms,f1,iPlt,titlestr)

if isstruct(res)
    RTs = res.response.RT';
    nTrials = length(res.response.RT);
    namestr = [res.params.initials ', ' num2str(res.params.tFreq) ' Hz Target'];
    bFit = true;
elseif isempty(res)
    bFit = false;
    nTrials = size(dict,1);
else
    RTs = res;
    nTrials = length(RTs);
    namestr = [];
    bFit = true;
end

%% Simulate RTs:
if size(weights,1)>size(weights,2)
    weights = weights';
end
w = repmat(weights,[nTrials 1]);
A = w.*dict;
simRTs = sum(A,2);

%% Find residuals & RMSE:
if bFit
    resid = (RTs-simRTs).^2;
    RMSE = sqrt(sum(resid)/length(RTs));
    
    
    R = corrcoef(RTs,simRTs);
    R = R(2);
    R2 = R^2;
    
    if bPlot
        if ~exist('f1','var') || isempty(f1)
            f1 = figure('color',[1 1 1],'pos',[0 0 600 600],'name',namestr);
        end
        if ~exist('iPlt','var'), iPlt = [1 1 1]; end
        
        subplot(iPlt(1),iPlt(2),iPlt(3)); hold on;
        plot(RTs,simRTs,'k.');
        
        % Plot unity line:
        maxRT = max([max(RTs) max(simRTs)]);
        minRT = min([min(RTs) min(simRTs)]);
        if minRT > 0, minRT = 0; end
        uu = minRT:maxRT;
        plot(uu,uu,'k');
        coef1 = polyfit(RTs,simRTs,1);
        set(gca,'FontSize',24);
        
        if exist('nTerms','var')
            title({titlestr;[num2str(nTerms) ' Terms (' num2str(nLinearTerms) ' L)'];...
                ['R^2 = ' num2str(R2,4)]},'FontSize',12);
        else
            title({['slope: ' num2str(coef1(1),4) '; int: ' num2str(coef1(2),4)];['R^2 = ' num2str(R2,4)]});
        end

        axis([minRT maxRT minRT maxRT]); axis square;
        xlabel('Observed RT (s)'); ylabel('Predicted RT (s)');
    else
        f1 = [];
    end
else
    RMSE=[];R2=[];f1=[];R=[];
end
end