% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL11_4vAcquireRF.m - Example of RF data acquisition
%
% Description:
%   Sequence programming file for L11-5v Linear array, acquiring RF data of
%   a single plane wave transmit and receive acquisition. All 128 transmit
%   and receive channels are active for each acquisition. External
%   processing is used asynchronous with respect to acquisition.
%
% Last update:
% 11/10/2015 - modified for SW 3.0
% 05/12/2020 - Update to SW 4.3 format for new user UIControls and External function definitions (VTS 1691)
%   (~/Example_Scripts/Vantage_Features/New UI Scheme/SetUpL11_5vFlash_NewUI)
clc;
clear all;

P.maxHighVoltage = 40;
P.sync=0;           % sw and hw are set to be synchronised if the value is set to 1
P.HV = 35;  

% Specify system parameters
Resource.Parameters.numTransmit = 256;      % no. of transmit channels (2 brds).
Resource.Parameters.numRcvChannels = 256;    % no. of receive channels (2 brds).
Resource.Parameters.speedOfSound = 1540;    % speed of sound in m/sec
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;       % runs script in simulate mode
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.
Resource.Parameters.Connector = [1 2];
Resource.Parameters.fakeScanhead = 1;

% Specify media points
Media.MP(1,:) = [0,0,100,1.0]; % [x, y, z, reflectivity]
Media.function = 'movePoints';

% Specify Trans structure array.
Trans.name = 'Doppler-3MHz-CR';
    Trans.units = 'mm'; 
    Trans.frequency = 3;
    Trans.Bandwidth = 3*[0.7, 1.3];% default assumed value of 60% of center frequency
    Trans.id = -1;
    Trans.type = 1;     % Array geometry is curved linear (x and z values only).
    Trans.connType = 1; % HDI connector
    Trans.numelements = 256;
    scanangle = 360 * (pi/180);    % degrees converted to radians
    radiusMm = 100;  % radius in mm.
    spacingMm = radiusMm * scanangle/(Trans.numelements); % spacing in mm.
    kerf = 0.2;   % guess (in mm)
        Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
        Trans.elevationFocusMm = 0; % nominal elevation focus depth from lens on face of transducer (unknown)
    Trans.radiusMm = radiusMm; % radius in mm.
    Trans.spacingMm = spacingMm;  % Spacing in mm.
    Trans.elementWidth = (spacingMm - kerf);  % width in mm
    deltatheta = scanangle/(Trans.numelements);
    firstangle = -pi + deltatheta/2; %   first element angle = -0.65405 radians
    lastangle = pi - deltatheta/2;
    %   Set default element positions (units in mm).
    Trans.ElementPos = zeros(Trans.numelements,5);
    Angle(1:256) = firstangle:deltatheta:lastangle;
    Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
    Trans.ElementPos(:,2) = 0;
    Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle);
    Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
    if ~isfield(Trans,'ElementSens')                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
        Theta = (-pi/2:pi/100:pi/2);
        Theta(51) = 0.0000001;                          % set to almost zero to avoid divide by zero.
        eleWidthWl = Trans.elementWidth * Trans.frequency/Resource.Parameters.speedOfSound;
        Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
    end
    Trans.impedance = 50; % using default value
    Trans.ConnectES = (1:256)';
    Trans.lensCorrection = 0;
    Trans.impedance = 50;
    Trans.maxHighVoltage = P.maxHighVoltage;
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000); % conversion factor from mm to wavelengths
% regardless of units, always provide spacing in wavelengths
Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths
Trans.radius = ceil(Trans.radiusMm * scaleToWvl);

% Specify Resource buffers.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 6400*256;  % this allows for 1/4 maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 1;       % allocate 10 frames.

% Now convert all units as required, based on Trans.units
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000); % conversion factor from mm to wavelengths
% regardless of units, always provide spacing in wavelengths
Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths

% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,10,1];


TX = repmat(struct('waveform', 1, ...
                   'Origin', zeros(1,3), ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,256), ...
                   'Delay', zeros(1,256)), 1, 256);

% - Set event specific TX attributes.
SA=1;
for n = 1:256-(SA-1)   % P.Elnum transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = scaleToWvl*Trans.ElementPos(n,1:3);
    % Set transmit Apodization so that only one element is active.
    TX(n).Apod(n:n+(SA-1)) = 1.0;    % Only one active transmitter for each TX.
    
end
for n = 256-(SA-1)+1:256   % P.Elnum transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = scaleToWvl*Trans.ElementPos(n,1:3);
    % Set transmit Apodization so that only one element is active.
    TX(n).Apod(n: 256) = 1.0;   
    TX(n).Apod(1:n-256+(SA-1)) = 1.0;    
end

TPC(1).hv = P.HV;
% Specify TGC Waveform structure.
TGC(1).CntrlPts = ones(1,8)*10;%[500,590,650,710,770,830,890,950];
TGC(1).rangeMax = 400;
TGC(1).Waveform = computeTGCWaveform(TGC);

A = [-0.00570 -0.00050 +0.01150 +0.01970 +0.00950 -0.02080 -0.04970 -0.04110 +0.02850 +0.14440 +0.25450 +0.29970];

B = [+0.00043 +0.00037 -0.00272 -0.00095 +0.00531 +0.00101 -0.00189 +0.00070 -0.01328 -0.00342 +0.03247 +0.00418 -0.02811 -0.00067 -0.02930 -0.00516 +0.13770 +0.00769 -0.24860 -0.00369 +0.29584];

% Specify Receive structure array -
Receive = repmat(struct(...
                'Apod', zeros(1,Trans.numelements), ...
                'startDepth', 0, ...
                'endDepth', TGC(1).rangeMax, ...
                'TGC', 1, ...
                'mode', 0, ...
                'bufnum', 1, ...
                'framenum', 1, ...
                'acqNum', 1, ...
                'LowPassCoef', A,...
                'InputFilter', B,...
                'sampleMode', 'NS200BW',...
                'callMediaFunc', 0), ...
                1,Resource.RcvBuffer(1).numFrames*256);
         

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive((i-1)*256+1).callMediaFunc = 1;
    for j = 1:256
    % -- full aperture acquisition.
        Receive((i-1)*256+j).Apod = [ones(1,128), ones(1,128)];%[ones(1,128), zeros(1,128)]
        Receive((i-1)*256+j).demodFrequency = TW(1).Parameters(1);
        Receive((i-1)*256+j).framenum = i;
        Receive((i-1)*256+j).acqNum = j;
    end
end

% Specify an external processing event.
Process(1).classname = 'External';
Process(1).method = 'myProcFunction';
Process(1).Parameters = {'srcbuffer','receive',... % name of buffer to process.
                         'srcbufnum',1,...
                         'srcframenum',-1,... % process the most recent frame.
                         'dstbuffer','none'};

% Specify sequence events.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;
SeqControl(2).command = 'jump';
SeqControl(2).argument = 1;
SeqControl(3).command = 'returnToMatlab';
nsc = 4; % nsc is count of SeqControl objects
lastTTHnsc = 0; % this variable keeps track of the last 'transferToHost' SeqControl index.
currentTTHnsc = 0; % this variable keeps track of the current 'transferToHost' SeqControl index.

n = 1;   % start index for Events
for i = 1:Resource.RcvBuffer(1).numFrames

    for j = 1:256                % Acquire frames
        Event(n).info = 'Acquisition.';
        Event(n).tx = j;
        Event(n).rcv = 256*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    Event(n-1).seqControl = nsc; % modify last acquisition Event's seqControl
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      if P.sync == 1
         SeqControl(nsc).condition = 'waitForProcessing';
         SeqControl(nsc).argument = lastTTHnsc;
         currentTTHnsc = nsc;
      end
      nsc = nsc+1;

	Event(n).info = 'Call external Processing function.';
	Event(n).tx = 0;
	Event(n).rcv = 0;
	Event(n).recon = 0;
	Event(n).process = 1;
	Event(n).seqControl = 3;
	if P.sync ==1
    Event(n).seqControl = [nsc,nsc+1,3];
       % The 'waitForTransferComplete' and 'markTransferProcessed' commands are
       % automatically executed with a reconstruction event, and don't need to be provided
       % by the user. They are duplicated here for helping user understanding the action
       SeqControl(nsc).command = 'waitForTransferComplete';
       SeqControl(nsc).argument = lastTTHnsc;
       nsc = nsc + 1;
       SeqControl(nsc).command = 'markTransferProcessed';
       SeqControl(nsc).argument = lastTTHnsc;
       nsc = nsc + 1;
       lastTTHnsc = currentTTHnsc;
    end
    n = n+1;
end

if P.sync ==1
    % fix the lastTTHnsc = 0 argument in SeqControl(4:6), change it to point to last TTH from acquisition loop
    SeqControl(4).argument = lastTTHnsc;
    SeqControl(5).argument = lastTTHnsc;
    SeqControl(6).argument = lastTTHnsc;
end

Event(n).info = 'Jump back to Event 1.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 2;

% User specified UI Control Elements and callback function for channel selection
nr = Resource.Parameters.numRcvChannels;
UI(1).Control = vsv.seq.uicontrol.VsSliderControl('LocationCode','UserB1','Label','Plot Channel',...
                  'SliderMinMaxVal',[1,256,128],...
                  'SliderStep',[1/nr,8/nr],'ValueFormat','%3.0f', ...
                  'Callback', @(~,~,UIValue)assignin('base','myPlotChnl',round(UIValue)) );

%External function
EF(1).Function = vsv.seq.function.ExFunctionDef('myProcFunction', @myProcFunction);

% Save all the structures to a .mat file.
save('MatFiles/RA_256_AcquireRF');

filename = 'MatFiles/RA_256_AcquireRF';
VSX;


%% **** Callback routines used by External function definition (EF) ****

function myProcFunction(RData)
    persistent myHandle
    Receive = evalin('base','Receive');
    % If myPlotChnl exists, read it for the channel to plot.
    if evalin('base','exist(''myPlotChnl'',''var'')')
        channel = evalin('base','myPlotChnl');
    else
        channel = 128;  % Channel no. to plot
    end
    % Create the figure if it does not exist.
    if isempty(myHandle)||~ishandle(myHandle)
        figure('name','Receive Signal','NumberTitle','off');
        myHandle = axes('XLim',[0,Receive(1).endSample],'YLim',[-2048 2048], ...
                        'NextPlot','replacechildren');
    end
    % Plot the element's RF data.
    V = version;
    MatlabV = V(end-5:end-1);
    plot(myHandle,RData(1:Receive(1).endSample,channel));
    if strcmp(MatlabV,'2014a') || str2double(MatlabV(1:end-1))<2014
        drawnow
    else
        drawnow limitrate
    end
    return
end
