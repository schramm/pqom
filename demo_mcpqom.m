%% read mocap data
d = mcread('PQoM_120BPM0002.mat');

%% mcpqom  parameters
markerNames1 = {'hand'}; % marker(s) to extract PQoM from
frameRange = []; % using full range
bpm = 120; % 120 beats per minute
noteDiv = '2 2t 4 8 8d 16' % note divisions (half, half triplets, quarter, eight, dotted eight, sixteenth)
windowSize = 10; % 10 beats
hopSize = .1; % 1/10 of a beat

%% extract PQoM and plot the results
[PQoM_Buffer, pqom_parameters] = mcpqom(d, markerNames1, frameRange, bpm, noteDiv, windowSize, hopSize); %compute PQoM
mcplotpqom(PQoM_Buffer, pqom_parameters); %plot PQoM
