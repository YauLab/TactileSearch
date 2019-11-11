function TopHammingOMP(res10,res30,res10ps,res30ps)

nBSIter = 1000;
sigLevel = 0.05;
sigLevel2 = 0.1;
ixCutoff = nBSIter*sigLevel;
ixCutoff2 = nBSIter*sigLevel2;
nSubj = size(res10,1);

ref = MakeTermReference;

% For each subject, concatenate aligned & scrambled data separately
% & compute HD between those 2 sequences:
obsHD = nan([1 nSubj]);
bshdAll = nan([nSubj nBSIter]);
pAll = false([1 nSubj]);
for ii = 1:nSubj
    
    t10 = PrepareTerms(res10(ii),ref);
    t30 = PrepareTerms(res30(ii),ref);
    t10ps = PrepareTerms(res10ps(ii),ref);
    t30ps = PrepareTerms(res30ps(ii),ref);    

    % Concatenate sync & phsc sequences together
    tsync = [t10 t30]; tphsc = [t10ps t30ps];
    % Find observed HD:
    obsHD(ii) = pdist([tsync;tphsc],'hamming');
    
    % Make distributions of permuted HDs & test:
    bsSeqSync = BSSequence2(t10,t30,nBSIter);
    [~,bshdAll(ii,:)] = BSSequence2(t10ps,t30ps,nBSIter,bsSeqSync);
    
    % Find p value of this obs HD within bs distribution:
    bs1 = sort(bshdAll(ii,:));
    cutoff = bs1(ixCutoff);
    pAll(ii) = obsHD(ii) < cutoff; % pAll indicates which subjects had significant obsHD
end

disp('Mean & sem of observed within-subject HDs:');
mhd1 = mean(obsHD)
sehd1 = std(obsHD)/sqrt(length(obsHD)) %#ok<NASGU>

% Compare observed mean HD across all comparisons to null meta
% distribution:
blank = nan([1 nSubj]);
nullMeanDist = nan([1 nBSIter]);
for jj = 1:nBSIter
    vals = blank;
    for kk = 1:nSubj
        vals(kk) = datasample(bshdAll(kk,:),1);
    end
    nullMeanDist(jj) = nanmean(vals(:));
end

% Determine where the observed average HD falls within the meta 
% null distribution:
bsVec = sort(nullMeanDist);
cutoff = bsVec(ixCutoff);
disp('Significantly different than chance?');
mhd1 < cutoff %#ok<*NOPRT>


%% Across subject comparison - concatenate OBSERVED & BS: 
bsAllSubj = nan([nSubj nSubj nBSIter]);
hdAllSubj = nan([nSubj nSubj]);
p1 = false([nSubj nSubj]); p2 = false([nSubj nSubj]);
bDone = false([nSubj nSubj]);
for ii = 1:nSubj
    t101 = PrepareTerms(res10(ii),ref);
    t301 = PrepareTerms(res30(ii),ref);
    t10ps1 = PrepareTerms(res10ps(ii),ref);
    t30ps1 = PrepareTerms(res30ps(ii),ref);    
    h1 = [t101 t301 t10ps1 t30ps1];  
    % Make bs sequences (permute then concatenate, x 1000):
    bsSeq1 = BSSequence(t101,t301,t10ps1,t30ps1,nBSIter);
    
    for jj = 1:nSubj
        if ~(bDone(ii,jj) || bDone(jj,ii))

            t102 = PrepareTerms(res10(jj),ref);
            t302 = PrepareTerms(res30(jj),ref);
            t10ps2 = PrepareTerms(res10ps(jj),ref);
            t30ps2 = PrepareTerms(res30ps(jj),ref);
            h2 = [t102 t302 t10ps2 t30ps2];
            
            % Find observed HD:
            hdAllSubj(ii,jj) = pdist([h1;h2],'hamming');
            
            % Make bs sequences (permute then concatenate, x 1000):
            [~,bsAllSubj(ii,jj,:)] = BSSequence(t102,t302,t10ps2,t30ps2,nBSIter,bsSeq1);
                        
            % Check where the observed HD falls in the null distribution:
            bs1 = sort(squeeze(bsAllSubj(ii,jj,:)));
            cutoff = bs1(ixCutoff);
            p1(ii,jj) = hdAllSubj(ii,jj) < cutoff;
            cutoff2 = bs1(ixCutoff2);
            p2(ii,jj) = hdAllSubj(ii,jj) < cutoff2;
            bDone(ii,jj) = true; bDone(jj,ii) = true;
        end
        if ii == jj
            hdAllSubj(ii,jj) = nan;
            p1(ii,jj) = false;
            p2(ii,jj) = false;
            bsAllSubj(ii,jj,:) = nan(size(bsAllSubj(ii,jj,:)));
        end
    end
end

hd1 = hdAllSubj(:); hd1 = hd1(~isnan(hd1));

disp('Mean & sem of observed across-subject comparison HDs:');
mhd1 = mean(hd1)
sehd1 = std(hd1)/sqrt(length(hd1)) %#ok<NASGU>

% Compare observed mean HD across all comparisons to null meta
% distribution:
blank = nan([nSubj nSubj]);
nullMeanDist = nan([1 nBSIter]);
for jj = 1:nBSIter
    bDone = false([nSubj nSubj]);
    vals = blank;
    for kk = 1:nSubj
        for ll = 1:nSubj
            if kk ~= ll && ~bDone(kk,ll)
                vals(kk,ll) = datasample(bsAllSubj(kk,ll,:),1);
            end
            bDone(kk,ll) = true; bDone(ll,kk) = true;
        end
    end
    nullMeanDist(jj) = nanmean(vals(:));
end

bsVec = sort(nullMeanDist);

% Determine where the observed average HD falls within the meta 
% null distribution:
cutoff = bsVec(ixCutoff);
disp('Significantly different than chance?');
mhd1 < cutoff %#ok<*NOPRT>

end

function tList = PrepareTerms(res,ref)
t = res.optBIC.terms;
if length(t) < 17, t = 1:17; end
% Turn term list into binary vector:
tList = false([1 length(res.terms)]);
tList(t) = true;
% Remove terms 1-17:
tList = double(tList(ref.NLTerms));
end

% Randomly permute 2 sequences & concatentate in nIter bs distribution:
function [bsSeq,hdAll] = BSSequence2(t1,t2,nIter,bsOther)
r = {t1,t2};
ll = [length(t1) length(t2)];
its = [0 cumsum(ll)];
bsSeq = nan([nIter its(end)]);

if exist('bsOther','var')
    bTest = true; 
    hdAll = nan([1 its(end)]); 
else
    bTest = false;
end

if bTest
    for ii = 1:nIter
        for kk = 1:length(r)
            bsSeq(ii,its(kk)+1:its(kk+1)) = shuffleVector(r{kk});
        end
        % Test hamming distance:
        hdAll(ii) = pdist([bsSeq(ii,:);bsOther(ii,:)],'hamming');
    end
else
    for ii = 1:nIter
        for kk = 1:length(r)
            bsSeq(ii,its(kk)+1:its(kk+1)) = shuffleVector(r{kk});
        end
    end
    hdAll = [];
end
end

% Randomly permute 4 sequences & concatentate in nIter bs distribution:
function [bsSeq,hdAll] = BSSequence(t10,t30,t10ps,t30ps,nIter,bsOther)
r = {t10,t30,t10ps,t30ps};
ll = [length(t10) length(t30) length(t10ps) length(t30ps)];
its = [0 cumsum(ll)];
bsSeq = nan([nIter its(end)]);

if exist('bsOther','var')
    bTest = true; 
    hdAll = nan([1 its(end)]); 
else
    bTest = false;
end

if bTest
    for ii = 1:nIter
        for kk = 1:length(r)
            bsSeq(ii,its(kk)+1:its(kk+1)) = shuffleVector(r{kk});
        end
        % Test hamming distance:
        hdAll(ii) = pdist([bsSeq(ii,:);bsOther(ii,:)],'hamming');
    end
else
    for ii = 1:nIter
        for kk = 1:length(r)
            bsSeq(ii,its(kk)+1:its(kk+1)) = shuffleVector(r{kk});
        end
    end
    hdAll = [];
end
end