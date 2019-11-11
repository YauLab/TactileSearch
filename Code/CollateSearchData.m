function resAll = CollateSearchData(initials,tFreq,dFreq,bPhJitter,homedir)

%% Initialize:
if ~exist('initials','var')
    initials = upper(input('Enter subject initials: '));
else 
    initials = upper(initials);
end
if ~exist('tFreq','var'), tFreq = 30; end
if ~exist('dFreq','var'), dFreq = 10; end
if ~exist('bPhJitter','var'), bPhJitter = false; end

%% Create filename:
datadir = fullfile(homedir,'Data',initials);
if bPhJitter
    basename = ['PBS_Jitter_t' num2str(tFreq) '_d' num2str(dFreq) '_s'];
else
    basename = ['PBS_t' num2str(tFreq) '_d' num2str(dFreq) '_s'];
end
sessNum = 1;
filename = [basename num2str(sessNum) '.mat'];
filepath = fullfile(datadir, filename);
if ~exist(filepath,'file')
    disp('No sessions collected yet! Returning...');
    resAll = [];
end

%% Collate & combine session results:
bFirst = true;
while exist(filepath,'file')
    % Load session results:
    load(filepath,'res');
    % Initialize w/ first session:
    if bFirst
        resAll = res;
        bFirst = false;
    else
        % Combine with previous session results:
        resAll = CombinePBS({resAll, res});
    end
    
    % Create filename of next session:
    sessNum = sessNum + 1;
    filename = [basename num2str(sessNum) '.mat'];
    filepath = fullfile(datadir, filename);
end
end