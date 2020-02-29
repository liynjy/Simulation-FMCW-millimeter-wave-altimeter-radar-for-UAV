%% Run fangdi_24G to generate simulated signals before running this program.

%% Load data from file, or new generate.
clear;clc
datafile = 'vrv.mat';  % select data file.

if ~exist('xrv')
    if exist(datafile)
        load(datafile);
    else
        fangdi_24G;
    end
end

%% Simulation of 2D-FFT processing on chirp.
% to 1024 points
dres = c/2/(bw*1024/1250);
xdata = xrv(101:1124,:);

% down sample
xd = xdata(1:4:end,:);

% range FFT
clear f fc
for k=1:32
    f = fft(xd(:,k).*hanning(256));
    fc(:,k) = f(1:32);
end

% doppler FFT -- chirp output data
for t=1:32
    fd(t,:)= fftshift(fft(fc(t,:).*hanning(32)'));
end

figure; 
subplot(221);mesh(abs(fd(1:10,:)));ylabel('range');xlabel('velocity');
subplot(222);imagesc(abs(fd(1:10,:)));ylabel('range');xlabel('velocity');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Algorithm
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% find max position index
absf = abs(fd);
[vidx,ridx] = max_pos_2D(absf);

%% RANGE
% focus vidx - range signal
rf = fd(:,vidx);

% back to fast time raw chirp
rf(64) = 0;
subplot(413);plot(abs(rf(1:32)))
xlabel('range');ylabel('Amplitude');xlim([0,32])

% add zeros, interpolation 
interpolation_N = 16;
iN = 2^nextpow2(interpolation_N);

dres = dres/iN;
ix = ifft(rf);
ix(64*iN) = 0;
irf = fft(ix);
subplot(414);plot(abs(irf(1:32*iN)))
xlabel('range');ylabel('Amplitude');xlim([0,32*iN])

% range
Range_Calib = dres*iN*0.5; % meters
R = max_pos_2D(abs(irf(1:32*iN))) - 1;
R = R*dres - Range_Calib;
disp(['R: ',num2str(R)])

%% VELOCITY
vres = lambda/2/(32*tm);
vf = fd(ridx,:);
V = center_gravity(abs(vf))- 17;
V = V*vres;
disp(['V: ',num2str(V)])

