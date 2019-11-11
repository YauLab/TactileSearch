function res = ModifyTargetIDs(res)

% Changed order of site IDs (in analyses)
%         from
% LE LI LK LF RF RK RI RE
%          to
% LE LI LK LF RE RI RK RF

if ~isempty(res)
    % Modify target IDs:
    b5 = res.params.iTarget == 5;
    b6 = res.params.iTarget == 6;
    b7 = res.params.iTarget == 7;
    b8 = res.params.iTarget == 8;
    res.params.iTarget(b5) = 8;
    res.params.iTarget(b6) = 7;
    res.params.iTarget(b7) = 6;
    res.params.iTarget(b8) = 5;
    
    % Modify distractor IDs:
    nTrials = length(res.params.iTarget);
    for ii = 1:nTrials
        ds = res.params.dSites{ii};
        b5 = ds == 5;
        b6 = ds == 6;
        b7 = ds == 7;
        b8 = ds == 8;
        ds(b5) = 8;
        ds(b6) = 7;
        ds(b7) = 6;
        ds(b8) = 5;
        res.params.dSites{ii} = ds;
    end
end
end