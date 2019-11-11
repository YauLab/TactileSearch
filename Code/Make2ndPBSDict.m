function [D, NLsites] = Make2ndPBSDict(res)
% D = Make2ndPBSDict(res)
%
% Make dictionary of all 1-site & 2-site terms.
%
% ejh 01/13/18

NLsites = combnk(2:17,2);
if ~exist('res','var') || isempty(res)
    D = []; return;
end

nTrials = length(res.params.iTarget);
nLTerms = 16; nNLTerms = size(NLsites,1);
nTerms = 1 + nLTerms + nNLTerms;
D = false(nTrials,nTerms);
D(:,1) = true;

% Fill dictionary, 1 trial at a time:
for ii = 1:nTrials
    % Find trial stimulus:
    iT = res.params.iTarget(ii);
    dSites = res.params.dSites{ii};
    
    % Deal with 1-target term:
    if iT > 0
        D(ii,iT+1) = true;
    end
    
    % Deal with 1-distractor terms:
    if ~any(dSites == 0)
        for dd = 1:length(dSites)
            D(ii,dSites(dd)+9) = true;        
        end
    end
end

% Deal with 2-site terms:
for jj = 1:size(NLsites,1)
    site1 = NLsites(jj,1); site2 = NLsites(jj,2);        
    D(:,jj+17) = D(:,site1)&D(:,site2);
end

end