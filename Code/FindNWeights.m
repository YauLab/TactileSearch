function [weights, D] = FindNWeights(res,terms,dict)

RT = res.response.RT';
% RT = res.response.RT(~isnan(res.response.RT));

% Remove all columns from dictionary that aren't being kept & reorder
% according to term order:
D = dict(:,terms);

% Normalize input RTs (divide by mean):
% meanY = nanmean(RT);
% y = RT'/meanY;

% Run linear regression:
% weights = regress(y,D);
weights = D\RT;

end