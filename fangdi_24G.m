% Author: Lin Junyang
% Date: 2019-11-09
% 1T1R Simulation. 
% Senario: UAV radar to horizontal ground/slope ground , height measurement.

clc;clear

%% Radar Parameters
fc = 24e9;   
c = physconst('LightSpeed');
lambda = c/fc;

tm = 5e-4;   % Chirp Cycle
bw = 300e6;    % FMCW Bandwidth
range_max = 5;     % Max detection Range 1~100 meters
v_max = 2.5;         % Max Velocity
%
range_res = c/2/bw;
sweep_slope = bw/tm;
fr_max = range2beat(range_max,sweep_slope,c);
fd_max = speed2dop(2*v_max,lambda);
fb_max = fr_max+fd_max;
fs = max(2*fb_max,bw);

%%
%% Use Phased Array System Toolbox to generate an FMCW waveform
waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,...
'SampleRate',fs);
%%
tx_antenna = phased.IsotropicAntennaElement('FrequencyRange',[23.8e9 24.4e9],'BackBaffled',true);
rx_antenna = phased.IsotropicAntennaElement('FrequencyRange',[23.8e9 24.4e9],'BackBaffled',true);
%%
transmitter = phased.Transmitter('PeakPower',0.001,'Gain',20);
receiver = phased.ReceiverPreamp('Gain',20,'NoiseFigure',8.5,'SampleRate',fs);

txradiator = phased.Radiator('Sensor',tx_antenna,'OperatingFrequency',fc,...
'PropagationSpeed',c);
rxcollector = phased.Collector('Sensor',rx_antenna,'OperatingFrequency',fc,...
'PropagationSpeed',c);

rng(2020);
fs_d = 2500000;
Dn = fix(fs/fs_d);

%%
%% --------------Radar Motion Platform-------------- %%
radar_s = phased.Platform('InitialPosition',[0;0;0],...
    'Velocity',[0.05;2.3;-0.04]);  %% *********** Set Radar Velocity Here **************        

%% Targets ------------- Ground -------------------- %%
target_ypos = -6:0.15:6;
target_num = size(target_ypos,2);
target_xpos = 1.3*ones(1,target_num) + 0*1.1*target_ypos; %% *********** Set Ground Shape Here ************** 
target_zpos = zeros(1,target_num);

target_pos = [[target_xpos,target_xpos,target_xpos];
    [target_ypos,target_ypos,target_ypos];
    [target_zpos-0.15,target_zpos,target_zpos+0.155]];
target_num = target_num*3;

target_rcs = 0.02*ones(1,target_num);
targets_vel = [zeros(1,target_num);zeros(1,target_num);zeros(1,target_num)];

targets = phased.RadarTarget('MeanRCS',target_rcs,'PropagationSpeed',c,'OperatingFrequency',fc);
targetmotion = phased.Platform('InitialPosition',target_pos,...
    'Velocity',targets_vel);

%%
%% Signal Propogation
% simulation of free space propagtion
channel = phased.FreeSpace('PropagationSpeed',c,...
'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%%
%%
% Generate Time Domain Waveforms of Chirps
% xr is the data received at rx array

Nsweep = 32;               % Number of Chirps (IF signal) of this simulation

chirp_len = fix(fs_d*waveform.SweepTime);
xr = complex(zeros(chirp_len,1,Nsweep));

disp('The simulation will take some time. Please wait...')
for m = 1:Nsweep
    if mod(m,1)==0
        disp([num2str(m),'/',num2str(Nsweep)])
    end
    
    % Update radar and target positions
    [radar_pos,radar_vel] = radar_s(waveform.SweepTime);
    [tgt_pos,tgt_vel] = targetmotion(waveform.SweepTime);
    [~,tgt_ang] = rangeangle(tgt_pos,radar_pos);
    
    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);
    
    % Toggle transmit element
    txsig = txradiator(txsig,tgt_ang);
    
    % Propagate the signal and reflect off the target
    txsig = channel(txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);
    txsig = targets(txsig);
    
    % Dechirp the received radar return
    rxsig = rxcollector(txsig,tgt_ang);
    rxsig = receiver(rxsig);
    dechirpsig = dechirp(rxsig,sig);
    
    % Decimate the return to reduce computation requirements
    for n = size(xr,2):-1:1
        xr(:,n,m) = decimate(dechirpsig(1:chirp_len*Dn,n),Dn,'FIR');
    end
end

range_res = range_res*size(dechirpsig,1)/Dn/size(xr,1);

%%
xrv = squeeze(xr);
save('vrv.mat',...
    'xrv','fc','fs_d','c','tm','bw','waveform','range_res',...
    'Nsweep','chirp_len','Dn','fb_max','lambda',...
    'v_max','range_max')

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Part II: Signal Processing                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if ~exist('xrv')
    load('vrv.mat');
end
    
% FFT points
nfft_r = 2^nextpow2(size(xrv,1));
nfft_d = 2^nextpow2(size(xrv,2));
nfft_mul = 2;
ra_res = range_res*size(xrv,1)/nfft_mul/nfft_r;

% RDM Algorithm
rngdop = phased.RangeDopplerResponse('PropagationSpeed',c,...
    'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs_d,...
    'RangeMethod','FFT','PRFSource','Property',...
    'RangeWindow','Hann','PRF',1/waveform.SweepTime,...
    'SweepSlope',waveform.SweepBandwidth/waveform.SweepTime,...
    'RangeFFTLengthSource','Property','RangeFFTLength',nfft_r*nfft_mul,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',nfft_d*nfft_mul,...
    'DopplerWindow','Hann');

% RD Map
[resp,r,sp] = rngdop(xrv);


% % Range-Energy Calibration
% for k=size(resp,1)/2+1:size(resp,1)
%     resp(k,:,:) = resp(k,:,:) * (k-size(resp,1)/2)^3;
% end

subplot(221);plotResponse(rngdop,squeeze(xrv));axis([-2*v_max 2*v_max 0 range_max-0.05])


%respmap = mag2db(abs(resp));
respmap = abs(resp);
respmap = avg_filter_2D(respmap,1);
subplot(222);mesh(respmap(nfft_r*nfft_mul/2+1:nfft_r*nfft_mul/2+1+30*nfft_mul,...
    :))
    %nfft_d*nfft_mul/2-12*nfft_mul:nfft_d*nfft_mul/2+12*nfft_mul))

subplot(413);plot(sum(respmap(nfft_r*nfft_mul/2+1:nfft_r*nfft_mul/2+1+30*nfft_mul,...
    nfft_d*nfft_mul/2-1:nfft_d*nfft_mul/2+2),2))










