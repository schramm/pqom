
%% mcpqom  settings
frameRange = []; % using full range
bpm = 60; % 60 beats per minute
noteDiv = '1 2 4' % note divisions (whole, half and quarter notes)
windowSize = 4; % 4 beats
hopSize = .5;   % 1/2 beat
% =======================================================


%% demo 1
% =======================================================
%% load MoCap structure (example 1);
load('example1.mat');  % load MoCap into 'data' variable
% define markers to extract PQoM
markerNames1 = {'52_bow_frog'}; % demo 1
%% extract PQoM 
[PQoM_Buffer, pqom_parameters] = mcpqom(data, markerNames1, [], bpm, noteDiv, windowSize, hopSize);
mcplotpqom(PQoM_Buffer, pqom_parameters);  
%% demo 2
% =======================================================
%% load MoCap structure 
load('example1.mat'); 
% define markers to extract PQoM
markerNames2 = {'08_torso_spine_(T5)', '52_bow_frog'}; % demo 2
noteDiv = [.5 1 2.2];
%% extract PQoM 
[PQoM_Buffer, pqom_parameters] = mcpqom(data, markerNames2, [], bpm, noteDiv, windowSize, hopSize);
mcplotpqom(PQoM_Buffer, pqom_parameters); 

%% demo 3
% =======================================================
%% load MoCap structure 
load('example1.mat'); 
% define markers to extract PQoM
markerNames3 = {'08_torso_spine_(T5)', '52_bow_frog', '05_L_torso_shoulder'}; % demo 3
noteDiv = '1d 2 4t' % note divisions (whole [dot], half and quarter [tripplet])
%% extract PQoM 
[PQoM_Buffer, pqom_parameters] = mcpqom(data, markerNames3, [], bpm, noteDiv, windowSize, hopSize);
mcplotpqom(PQoM_Buffer, pqom_parameters); 



