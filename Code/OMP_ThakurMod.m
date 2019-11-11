function [xAll,bAll] = OMP_ThakurMod(D,RT,nTerms,keepTerms,nIter,w17)
% [xAll,bAll] = OMP_ThakurMod(D,RT,nTerms,keepTerms,nIter,w17)

if ~exist('nIter','var')
    nIter = 25;
end
cutThresh = 0.75;

[~,K] = size(D);

% Pre-allocate:
xAll = nan([nTerms nIter]);
bAll = nan([nTerms nIter]);

for iter = 1:nIter
    
    r = RT;              
    b = zeros(nTerms,1);        
    D_col = [];          
    
    % Determine if you're starting with selected terms already:
    if ~exist('keepTerms','var') || isempty(keepTerms)
        cnt = 0;
    elseif exist('w17','var')
        % Add 16 linear coefficients & constant to coefficient list:
        b = 1:17;
        % Find contribution of all 17 terms to this signal:
        D17 = D(:,1:17);
        x17 = D17 \ RT;
        % Find residual signal:   --> don't include 1:17 in future coefs
        r = RT - D17 * x17;
        % Update count:
        cnt = length(keepTerms);
    else
        cnt = length(keepTerms);
        bTmp = zeros(cnt,1);
        % Determine remaining signal (after inclusion of these terms):
        for ii = 1:cnt
            inds = setdiff(1:cnt,bTmp);
            x_tmp = zeros(cnt,1);
            for jj = inds
                x_tmp(jj) = D(:,jj)' * r /norm(D(:,jj));
            end
            [~,ichosen] = max(abs(x_tmp));
            bTmp(ii) = keepTerms(ichosen);
            D_col = [D_col D(:,bTmp(ii))];
            x_ls = D_col \ RT;
            r = RT - D_col * x_ls; % update residual
        end
        b(1:cnt) = bTmp;
    end
    
    while (cnt < nTerms)  
        cnt = cnt+1;
        x_tmp = nan(K,1);
        inds = setdiff(1:K,b);
        for i = inds
            x_tmp(i) = D(:,i)' * r / norm(D(:,i));
        end
        
        % At each iteration, randomly choose a term that has >75% of max proj:
        absP = abs(x_tmp);
        [maxP,ichosen] = nanmax(absP);
        cutoff = cutThresh*maxP;
        possIxs = find(absP >= cutoff);
        if ~isempty(possIxs)
            ichosen = datasample(possIxs,1);
        end
        
        b(cnt) = ichosen;
        D_col = [D_col D(:,ichosen)];                        %#ok<*AGROW>
        x_ls = D_col \ RT;     
        if exist('w17','var')
            % Add 17 fixed coefficients back in:
            D2 = [D17 D_col];
            x2 = [x17; x_ls];
            % Update residual:
            r = RT - D2 * x2;
        else
            % Update residual:
            r = RT - D_col * x_ls;       
        end
    end
    
    if exist('w17','var'), x_ls = x2; end
    xAll(:,iter) = x_ls;
    bAll(:,iter) = b;
end

end