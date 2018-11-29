function [pqom_buffer, pqom_parameters] = mcpqom(mocapData, markerName, frameRange, bpm, noteDiv, winSize, hopSize)
% Estimates the periodic quantity of motion from a MoCap structure.
%
% syntax
% [pqom_buffer, pqom_parameters] = mcpqom(mocapData, markerName, frameRange, bpm, noteDiv, winSize, hopSize)
% 
% input parameters
% mocapData: MoCap data structure.
% markerName: cell array containg the name (string) of skeleton markers where the PQoM will be estimated. 
%            If only one marker name is passed then the anchor reference is the origin [0 0 0]' of the world coordinate system. 
%            In case of two or more marker names are passed then the first marker name is the anchor (world reference) point.
%            eg.: markerName = {'torso_spine', 'left_hand'}; % The PQoM is estimated from 
%            the 'left_hand' marker using the 'torso_spine' as anchor point.
% frameRange: [ini end] array that specifies the range of data frames from mocapData (default: full range == [])
% bpm: beats per minute of the analised song (default: 60)
% noteDiv: string with the note divisions used on the PQoM estimation.
%            [1: whole note; 2: half note; 4: quarter note; 8: eight note; 16; 32; ] 
%            (default: '1 2 4 8')
%            Note divisions supports triplets and dotted notes.  Add 't' after the
%            number for triplets or 'd' for dotted notes (eg. '1 2d 4t ').
%            note divisions can be specified as a numeric array. In this case, 
%            the array values are interpreted as frequency in Hz. 
%            (eg. [3.5  4  5.2] will estimate PQoM on 3.5Hz, 4Hz and 5.2Hz)
% winSize: window size used to evaluate the PQoM (in beat times, default:  4*bpm)
% hopSize: the step  size between the windows of evaluation (in beat times, default:  1*bpm)  
%                    
% output
% pqom_buffer: buffer containg the estimated PQoM for each frequency band
% pqom_parameters: matlab structure containing the configuration parameters, 
%            including an array with the center frequency of each filter band used to estimate the PQoM 
%   
% example
% see demo_mcpqom.m
%
% comments
% The output variables might be used as input parameters to generate
%    plots of the PQoM estimations using the function 'mcplotpqom'. 
%         Ex.: mcplotpqom(pqom_buffer, pqom_parameters);  
%    If no ouputs are defined then the mcpqom generates plots by default.
%
% see also
% mcplotpqom
%
% ===============================================
% Version 1.0
% July 2015
% ===============================================
% Rodrigo Schramm (rodrigo.schramm@gmail.com)  and   Federico Visi (federico.visi@plymouth.ac.uk) 
% UFRGS / Brazil          ICCMR / Plymouth University / UK
% 
% If you plan to use this script, please cite the paper: 
%     Visi, Federico and Schramm, Rodrigo and Miranda, Eduardo. Gesture in Performance with Traditional Musical Instruments and Electronics: 
%     Use of Embodied Music Cognition and Multimodal Motion Capture to Design Gestural Mapping Strategies.
%     Proceedings of the 2014 International Workshop on Movement and
%     Computing. MOCO '14, p. 100-105, ACM, Paris, 2014.
% ===============================================
%
% Part of the Motion Capture Toolbox, Copyright 2008, 
% University of Jyvaskyla, Finland



%%% alloc data structures
%===============================================
qom = [];
fWn =[];
pqom_buffer = [];
%===============================================

%%% checking input arguments
%===============================================
if nargin<7
    hopSize = 1; 
end
if nargin<6
    winSize = 4;
end
if nargin<5
    noteDiv = [1 2 4 8]; 
end
%noteDiv = unique(noteDiv); % duplicated entries do not make sense
if nargin<4
    bpm = 60;
end
if nargin<3
   frameIni = 1;
   frameEnd = size(mocapData.data,1);
else 
    if isempty(frameRange) || length(frameRange)~=2
        frameRange = [1 size(mocapData.data,1)];
    end
    frameIni = frameRange(1);
    frameEnd = frameRange(2);
end
if nargin<2
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    error([10, 'Input arguments "mocapData" and "markerName" should be defined.', 10]);
end
%===============================================

%%% check frame range
%===============================================
if frameIni<1   
    frameIni=1;
    warning('mcpqom:input_parameters', 'Initial data frame is out of range. Using: frameIni=%d', frameIni);
elseif frameIni>size(mocapData.data,1) 
    frameIni = size(mocapData.data,1)-1;
    warning('mcpqom:input_parameters', 'Initial data frame is out of range. Using: frameIni=%d', frameIni);
end
if frameEnd<1   
    frameEnd=1;
    warning('mcpqom:input_parameters', 'Last data frame is out of range. Using: frameEnd=%d', frameEnd);
elseif frameEnd>size(mocapData.data,1)    
    frameEnd = size(mocapData.data,1);
    warning('mcpqom:input_parameters', 'Last data frame is out of range. Using: frameEnd=%d', frameEnd);
end

%===============================================

%%% extracting data range
%===============================================
data = mocapData.data(frameIni:frameEnd,:);
%===============================================

%%% filter settings
%===============================================
fs = mocapData.freq;
%ws = ceil( (60/bpm) * winSize * fs); % window size (in frames)
ws = max(round((60/bpm)*winSize*fs),1); % window size (in frames)
hs = max(round((60/bpm)*hopSize*fs),1); % hop size (in frames)

if ws>=size(data,1)
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    error('mcpqom:input_parameters', 'Input parameter [frameRange] is too small for the specified winSize [winSize=%2.2f]. \nframeRange size should have at least %d frames.',winSize, ws+1);    
end

anchor = [];
c = 1;
%%% iterate and get PQoM for each skeleton marker in markerName array.
for k=1:length(markerName)    
    mName = markerName{k}; % input parameter array
    idx = strcmp(mName, mocapData.markerName);
    if sum(idx)==0
        [y,fs] = audioread('mcsound.wav');
        sound(y,fs);
        error('mcpqom:input_parameters', 'MarkerName [%s] is not available in MoCap structure.', mName);
        return;
    end
    
    %%% get marker index
    idx =find(idx);
    idx = idx(1);
    
    % get the x,y,z movement 
    x = data(:,idx*3-2);
    y = data(:,idx*3-1);
    z = data(:,idx*3-0);

    % remove any NaN (//TODO)
    if sum(isnan(x))>0 | sum(isnan(y))>0 | sum(isnan(z))>0 
       warning('mcpqom:input_parameters','Part of signal contains NaN elements'); 
    end
    x(isnan(x))=0;
    y(isnan(y))=0;
    z(isnan(z))=0;

    %%% 
    if c==1 && length(markerName)>1  
        anchor = [x,y,z];
    elseif length(markerName)==1  %% anchor point is the origin [0 0 0]'.
        dd = ([x,y,z]).^2;
        dd = sqrt(sum(dd'))';
        qom(:,1) =  dd;        
        warning('mcpqom:input_parameters','Anchor point is [0 0 0].');
    else
        dd = (anchor - [x,y,z]).^2;
        dd = sqrt(sum(dd'))';
        qom(:,c-1) =  dd;
    end    
    c=c+1;
end


% padding signal to perform convolution
%===============================================
qom = [qom; qom(end:-1:end-ws,:)]; 

%%% filtering the signal according to specific note divisions (freq bands)
%===============================================
[pqom_buffer, fWn, fHz noteType] = filterBands(qom, fs, bpm, noteDiv, ws, hs);

% define output parameters for use in mcplotpqom function
pqom_parameters =[];
pqom_parameters.fWn=fWn;
pqom_parameters.fHz=fHz;
pqom_parameters.noteType = noteType;
pqom_parameters.bpm = bpm;
pqom_parameters.fs = fs;
pqom_parameters.noteDiv = noteDiv;
pqom_parameters.ws = ws;
pqom_parameters.hs = hs;
pqom_parameters.markerName = markerName;

if nargout==0
    mcplotpqom(pqom_buffer, pqom_parameters); 
end

end



function [buffer, freqWn, freqHz, noteType] = filterBands(skelData, fs, bpm, noteDiv, ws, hs)
%%% freq bands (based on noteDiv)
% [1: whole note; 2: half note; 4: quarter note; 8: eight note; 16; 32; ] 
% WARNING: 1 beat time == quarter note 
%===============================================
[fB, noteType] = convertNoteDiv2FreqHz(noteDiv,bpm);

% alloc memory
skelFilt = zeros([size(skelData) length(fB)]);
s = size(skelFilt(:,1,1));
buffer = zeros([length(1:hs:s(1)), length(fB)]);
freqWn = []; % cuttof frequency per band


%%% converts note duration to frequency (Hz)
for k=1:length(fB)
    freqWn(k) = 2*fB(k)/fs; % cutoff freq
    Wn = 2*[fB(k)-(fB(k)/10) fB(k)+(fB(k)/10)]/fs; % band pass 

    for j=1:size(skelData,2); % each marker
        v = skelData(:,j); 
        %sl = length(v);    
        %%% butterworth filter (second order, band pass)   
        [bf, af] = butter(2, Wn);
        out = filtfilt(bf,af,v); 
        out = out.^2;
        out = out/sum(out);
        skelFilt(:,j,k) = out;    
    end
end

%%% integrate along the window and markers, for each band
for k=1:length(fB)
    c=1;    
    for i=1:hs:s(1)-ws;
        v = sum(skelFilt(:,:,k),2); % sum along markers
        v = v(i:i+ws);
        ss = sum(v);  % sum along window
        buffer(c, k) = ss;        
        c=c+1;    
    end
end
              
%%% trim buffer tail (remove the padding)
buffer = buffer(1:end-ws/hs, :);
freqHz = freqWn*fs/2;
end


function [freqHz, noteType] = convertNoteDiv2FreqHz(noteDiv, bpm)
noteType = {};
freqHz = [];
if ischar(noteDiv)
    freqHz=[];
   r = noteDiv;
   k=0;
   while ~isempty(r)
       k=k+1;
        [n,r] = strtok(r);        
        tripleN = 1;
        dotN = 1;
        
        [b,t] = strtok(n, 't');
        if ~isempty(t)
            n = b;           
            tripleN = 2/3;
            t = '(t)';
        end
        [b,d] = strtok(n, 'd');
        if ~isempty(d)
            n = b;
            dotN = 1.5;  
            d = '(d)';
        end
        % //TODO implement case: 4T. ???
       n = str2num(n); 
        
    switch n
        case 1 % whole note;
            dur = 4;
            noteType{k}='whole note';
        case 2 % half note;
            dur = 2;
            noteType{k}='half note';
        case 4 % quarter note;
            dur = 1;
            noteType{k}='quarter note';
        case 8 % eight note;
            dur = 1/2;
            noteType{k}='eight note';
        case 16 % 16th note;
            dur = 1/4;
            noteType{k}='16th note';
        case 32 % 32th note;
            dur = 1/8;
            noteType{k}='32th note';
        otherwise 
            [y,fs] = audioread('mcsound.wav');
            sound(y,fs);
            error('mcpqom:input_parameters', 'Invalid note division. noteDiv==%s  \nThe string format only suport [1: whole note; 2: half note; 4: quarter note; 8: eight note; 16 or 32; ]', num2str(n));
    end
    %[dur, tripleN, dotN, dur*tripleN*dotN, 1/(dur*tripleN*dotN)],
    dur = dur*tripleN*dotN;
      noteType{k} = [noteType{k}, '', t, '', d,'']; 
      freqHz(k) = 1/(dur*(60/bpm));
   end
else % when noteDiv means Frequency (in Hz)
    freqHz = noteDiv;
    for k=1:length(noteDiv)
        noteType{k} = [num2str(noteDiv(k)), 'Hz'];
    end
end
end
