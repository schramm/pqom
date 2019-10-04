function [pqom_buffer, pqom_parameters] = mcpqom(mocapData, markerName, frameRange, bpm, noteDiv, winSize, hopSize)
% Estimates the periodic quantity of motion from a MoCap structure.
%
% syntax
% [pqom_buffer, pqom_parameters] = mcpqom(mocapData, markerName, frameRange, bpm, noteDiv, winSize, hopSize)
% 
% input parameters
% mocapData: MoCap data structure. 
% markerName: cell array containg the name (string) of skeleton markers where the PQoM will be estimated. 
%            In case of two or more marker names are passed then the PQoM will be the sum of PQoM estimates for each marker.
%            eg.: markerName = {'torso_spine', 'left_hand'}; 
%% TODO allow the use of both numbers and strings.
%
%%
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
%
% Updates: 
% Version 1.1 (in progress)
% August 2019

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

[freqCenterHz, noteType] = convertNoteDiv2FreqHz(noteDiv, bpm);

%%% settings
%===============================================
fs = mocapData.freq;
ws = max(round((60/bpm)*winSize*fs),1); % window size (in frames)
hs = max(round((60/bpm)*hopSize*fs),1); % hop size (in frames)


if ws>=size(data,1)
    [y,fs] = audioread('mcsound.wav');
    sound(y,fs);
    error('mcpqom:input_parameters', 'Input parameter [frameRange] is too small for the specified winSize [winSize=%2.2f]. \nframeRange size should have at least %d frames.',winSize, ws+1);    
end

c=0;
%%% iterate and get PQoM for each skeleton marker in markerName array.
for k=1:length(markerName)    
    mName = markerName{k}; % input parameter array
    idx = strcmp(mName, mocapData.markerName);
    if sum(idx)==0        
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

    % estimate PQoM at each axis (x,y and z)
    [q_x,ffq_x, ffHz_x,pqom_x] = rs_pqom(x', ws, hs, fs, freqCenterHz);
    [q_y,ffq_y, ffHz_y,pqom_y] = rs_pqom(y', ws, hs, fs, freqCenterHz);
    [q_z,ffq_z, ffHz_z,pqom_z] = rs_pqom(z', ws, hs, fs, freqCenterHz);       
    
    if (isempty(pqom_buffer)) 
        pqom_buffer = zeros(size(pqom_x));
        qom = zeros(size(q_x));
    end
   pqom_buffer =  pqom_buffer + pqom_x + pqom_y +pqom_z;        
   qom = qom + q_x + q_y + q_z;
   c = c+3; 
end

pqom_buffer = pqom_buffer/c;
qom = qom/c;

noteType{length(noteType)+1} = "residual";
   
% define output parameters for use in mcplotpqom function
pqom_parameters =[];
pqom_parameters.fWn=fWn;
pqom_parameters.fHz=[freqCenterHz, Inf];
pqom_parameters.noteType = noteType;
pqom_parameters.bpm = bpm;
pqom_parameters.fs = fs;
pqom_parameters.noteDiv = noteDiv;
pqom_parameters.ws = ws;
pqom_parameters.hs = hs;
pqom_parameters.markerName = markerName;
pqom_parameters.qom = qom;
pqom_parameters.ffHz_x = ffHz_x;
pqom_parameters.ffHz_y = ffHz_y;
pqom_parameters.ffHz_z = ffHz_z;

if nargout==0
    mcplotpqom(pqom_buffer, pqom_parameters); 
end

end

function [q,ffq, ffHz,pqom] = rs_pqom(sig, ws, hs, Fs, freqCenterHz)

[q,ffq,ffHz] = rs_qom(sig, ws, hs, Fs, freqCenterHz);
qq = repmat(q',[1,size(ffHz,2)]);
pqom = ffHz.*qq;

end


%% compute the QoM and the FFT analysis (per frame)
function [q,ffq, ffHz] = rs_qom(sig, winSize, hopSize,Fs, freqCenterHz)

if winSize>length(sig)
    error('signal length must be bigger than winSize');
end

s = length(sig);
hWinSize = floor(winSize/2);
% append signal
asig = [zeros([1 hWinSize])+sig(1) sig zeros([1 hWinSize])+sig(end)];
% alloc ouput
q = zeros([1 ceil(length(sig)/hopSize)]);

nPointFFT = max( 2^10, winSize); %% TODO

ffq = zeros([ceil(length(sig)/hopSize), nPointFFT/2+1]);
ffHz = zeros([ceil(length(sig)/hopSize), length(freqCenterHz)+1]);
c=0;
for i=hWinSize+1:hopSize:s+hWinSize
    ss = asig(i-hWinSize:i+hWinSize); % crop window
    d = diff(ss);                       % diff
    d = sum(abs(d));                    % TODO: if winSize is even, window size is equal to winSize+1
    c=c+1;
    q(c) = d / (winSize/Fs);  % QoM
    %q(c) = d;  % QoM
    
    %% FFT    
    %f = Fs*(0:(L/2))/L;    
    %n = round(freqHz*L/Fs);        
    [X] = rs_fft(ss,nPointFFT);        
    ffq(c,:) = X;
    
    [pks,lcs]=findpeaks(X, 'npeaks', 10); %//TODO how to specify the number of peaks?         
    n = round(freqCenterHz*nPointFFT/Fs); % center frequencies
    for j=1:length(n)        
        ffHz(c,j) = X( round(n(j)) );                         
    end
        
    % compute residual (sum of peaks that not belong to the freq bands we are analysing
    toDelete = zeros(size(lcs));
    for j=1:length(pks)        
        for k=1:length(n)
            if (abs(lcs(j)-n(k))<5)
                toDelete(j)=1;
            end
        end
    end
    pks(toDelete==1)=[];    
    ffHz(c,end) = sum(pks);   %% residual
    ffHz(c,:) = ffHz(c,:)/sum(ffHz(c,:)); % normalise
end


end

%% compute the FFT
function [X] = rs_fft(x, nPointFFT)
L = length(x);
h = hann(L)';
X = fft(x.*h, nPointFFT); % force fft with higher resolution (fft does pad with zeros)
X = abs(X/L);
P1 = X(1:nPointFFT/2+1);
%P1(2:end-1) = 2*P1(2:end-1); 
X = P1;
end

%% convert the periodic note duration/division (rhythm or beat) to Hz
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
%            [y,fs] = audioread('mcsound.wav');
%            sound(y,fs);
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

    