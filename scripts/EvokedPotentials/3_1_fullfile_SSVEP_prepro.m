%written by Renate Kat
%modified by FGO (5-2-20)

clear all
close all

[file,path]=uigetfile;

%load data
[header, signal_header, signal_cells] = blockEdfLoad(fullfile(path,file));

% Call the annotation conversion function to convert annotations to a list
% of event times. 
% event_times is an array of times recorded for rising or falling edges (in seconds)
% event_values are the integer (combined value) for the synchronisation 12-bit state
% event_binary are the binary (per input value) for the synchronisation for each of the 12-bits

[ event_times, event_values, event_binary ] = parse_taini_edf_annotations( signal_header, signal_cells );

fs = signal_header(1).samples_in_record / header.data_record_duration;
%%
%Convert event_times from seconds to datapoint by multiplying with the
%sampling frequency. Select every second points to select the 3's.
%transpose the array and remove first and last trigger. Round to end up at
%specific measurement points.
rSfreq=round(fs);

event_times([1]) = []; % removes the first trigger which usually marks the beginning of the recording
event_points = event_times * rSfreq; 
event_points = event_points(1, 1:2:end)';
event_points = round(event_points);
plot(event_points, 1, 'r*')
hold on
ylim([0.5 1.5])
plot(event_points, 0.99, 'b*')
for ii = 1:length(event_points)
    text(event_points(ii),0.97,num2str(ii),'Color','k')
end

%% Divide over right frequencies
trigV1 = event_points(1426:1475,:); % 2; only in case datapoint 1 causes a problem
trigV2 = event_points(1476:1525,:);
trigV3 = event_points(1526:1575,:);

%select the EEG right channels from the struct and downsample to rSfreq
Ch3 = cell2mat(signal_cells(1,3))';
Ch4 = cell2mat(signal_cells(1,4))';
Ch12 = cell2mat(signal_cells(1,12))';
Ch13 = cell2mat(signal_cells(1,13))';
Ch14 = cell2mat(signal_cells(1,14))';
Ch8 = cell2mat(signal_cells(1,8))';
Ch11 = cell2mat(signal_cells(1,11))';

EEG = [Ch3;Ch4;Ch12;Ch13;Ch14;Ch8;Ch11];

figure(2)
plot(EEG(1,:))

% Cut around triggers (norm -4000:4000)
stimMatV1 = bsxfun(@plus,trigV1,-4000:4000); 
stimMatV2 = bsxfun(@plus,trigV2,-4000:4000);
stimMatV3 = bsxfun(@plus,trigV3,-4000:4000);

stimMat = cat(3,stimMatV1', stimMatV2', stimMatV3');

T = (-4000*1/rSfreq):1/rSfreq:(4000*1/rSfreq);

%% Filtering and artifact rejection
for iCh = 1:size(EEG,1)

    %filter EEG data
    order = 4;
    hpFilt  = 0.5;
    lpFilt = 70;
    
    
    %Get rid of packeloss artificial low values and interpolate
    EEGCh = EEG(iCh,:);
    EEGCh (EEGCh < 6) = nan;
    EEGCh (EEGCh > 10) = nan;
    Packloss(1,iCh) = (sum(isnan(EEGCh))/size(EEGCh,2)) *100;
    EEGres = resample(EEGCh,1:length(EEGCh));
    
    %Filter EEG data
    EEGHp = filterEEG(EEGres, order, hpFilt, rSfreq, 'high', 2);
    EEGHLp = filterEEG(EEGHp, order, lpFilt, rSfreq, 'low', 2);
    
            
    %reject artifacts
    average = mean(EEGHLp);
    stdev = std(EEGHLp);
    EEGHLp(EEGHLp < average - 3*stdev) = nan;
    EEGHLp(EEGHLp > average + 3*stdev) = nan;
    Packloss(2,iCh) = (sum(isnan(EEGHLp))/size(EEGHLp,2)) *100;
    EEGHLp = resample(EEGHLp,1:length(EEGHLp));
    
    %Cut around triggers, average and baseline correct
   for iVolt = 1: size(stimMat,3)
       erp =(EEGHLp(stimMat(:,:,iVolt)));
       baseline = mean(erp(2000:3000));
       ERP(:,:,iVolt,iCh) = erp - baseline; 
    end
    
    
end

%% CONVOLUTION

% frequency parameters
min_freq =  5; % in Hz
max_freq = 70; % in HZ
num_freq = 65; % in count

% set range for variable number of wavelet cycles
range_cycles = [ 1 25 ];

freq = linspace(min_freq,max_freq,num_freq);
nCycs = linspace(range_cycles(1),range_cycles(end),num_freq);

% create a time window for complex Morlet wavelet
time = (0:2*rSfreq)/rSfreq;
time = time - mean(time); 

%put all trials after each other
ERPV1R = reshape(ERP(:,:,1,:),1,[]);
ERPV2R = reshape(ERP(:,:,2,:),1,[]);
ERPV3R = reshape(ERP(:,:,3,:),1,[]);

%% First stimulus paradigm
% Step 1: N's of convolution
ndata = length(ERPV1R);% note the different variable name!
nkern = length(time);
nConv = ndata + nkern -1;
halfK = floor(nkern/2);

%Fourier of the data (outside the loop cause the same every time)
dataXV1 = fft( ERPV1R ,nConv );
% initialize TF matrix
tfV1 = zeros(num_freq,length(T),size(ERP,4));
for fi=1:num_freq
    % create wavelet
    s = nCycs(fi)/(2*pi*freq(fi));
    
    % cmw = complex Morlet wavelet
    cmw  = exp(1i*2*pi*freq(fi).*time) .* exp(-time.^2./(2*s^2));
    % FFT
    cmwX = fft(cmw,nConv);
    % Normalizing
    cmwX = cmwX./max(cmwX);
    
    % the rest of convolution
    as = ifft( dataXV1.*cmwX );
    as = as(halfK+1:end-halfK);
    as = reshape(as,size(ERP(:,:,1,:)));                 
    
    % extract power
    aspow = abs(as).^2;
    % average over trials and put in matrix
    tfV1(fi,:,:) = mean(squeeze(aspow),2);                            
    
    tfV1(fi,:,:) = 10*log10( bsxfun(@rdivide, squeeze(tfV1(fi,:,:)), mean(squeeze(tfV1(fi,2000:3000,:)),1)) ); % interval of baseline
    
end

figure(8), clf
subplot(3,1,1)
contourf(T,freq,tfV1(:,:,1),40,'linecolor','none')
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Primary Visual')
subplot(3,1,2)
contourf(T,freq,tfV1(:,:,2),40,'linecolor','none')
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Primary Auditory')
subplot(3,1,3)
contourf(T,freq,tfV1(:,:,3),40,'linecolor','none')
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Prefrontal Cortex')
suptitle('10Hz')

%% second stimulus paradigm
% Step 2: N's of convolution
%ndata2 = length(ERPV2R);% note the different variable name!
%nkern = length(time);
%nConv2 = ndata2 + nkern -1;
%halfK = floor(nkern/2);

dataXV2 = fft( ERPV2R ,nConv );
tfV2 = zeros(num_freq,length(T),size(ERP,4));

for fi=1:num_freq
    
    % create wavelet
    s = nCycs(fi)/(2*pi*freq(fi));
    
    cmw  = exp(1i*2*pi*freq(fi)*time) .* exp(-time.^2./(2*s^2));
    
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % the rest of convolution
    as = ifft( dataXV2.*cmwX );
    as = as(halfK+1:end-halfK);
    as = reshape(as,size(ERP(:,:,1,:)));    
    % extract power
    aspow = abs(as).^2;
    
    % average over trials and put in matrix
    tfV2(fi,:,:) = mean(aspow,2);
    
    tfV2(fi,:,:) = 10*log10( bsxfun(@rdivide, squeeze(tfV2(fi,:,:)), mean(squeeze(tfV2(fi,2000:3000,:)),1)) );
    
end

figure(9), clf
subplot(3,1,1)
contourf(T,freq,tfV2(:,:,1),40,'linecolor','none')
%set(gca, 'clim', [-3,7.5],'xlim',[-2,5]) 
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Primary Visual')
subplot(3,1,2)
contourf(T,freq,tfV2(:,:,2),40,'linecolor','none')
%set(gca,'clim', [-3,7.5],'xlim',[-2,5]) %'clim', [-3,7.5]
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Primary Auditory')
subplot(3,1,3)
contourf(T,freq,tfV2(:,:,3),40,'linecolor','none')
%set(gca,'clim', [-3,7.5],'xlim',[-2,5]) %
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Prefrontal Cortex')
suptitle('20Hz')

%% third stimulus paradigm
% Step 2: N's of convolution
%ndata3 = length(ERPV3R);% note the different variable name!
%nkern = length(time);
%nConv3 = ndata3 + nkern -1;
%halfK = floor(nkern/2);

dataXV3 = fft( ERPV3R ,nConv );
tfV3 = zeros(num_freq,length(T),size(ERP,4));

for fi=1:num_freq
    
    % create wavelet
    s = nCycs(fi)/(2*pi*freq(fi));
    
    cmw  = exp(1i*2*pi*freq(fi)*time) .* exp(-time.^2./(2*s^2));
    
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(cmwX);
    
    % the rest of convolution
    as = ifft( dataXV3.*cmwX );
    as = as(halfK+1:end-halfK);
    as = reshape(as,size(ERP(:,:,1,:)));    
    % extract power
    aspow = abs(as).^2;
    
    % average over trials and put in matrix
    tfV3(fi,:,:) = mean(aspow,2);
    
    tfV3(fi,:,:) = 10*log10( bsxfun(@rdivide, squeeze(tfV3(fi,:,:)), mean(squeeze(tfV3(fi,2000:3000,:)),1)) );
end

figure(10), clf
subplot(3,1,1)
contourf(T,freq,tfV3(:,:,1),40,'linecolor','none')
%set(gca, 'clim', [-3,7.5],'xlim',[-2,5]) 
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Primary Visual')
subplot(3,1,2)
contourf(T,freq,tfV3(:,:,2),40,'linecolor','none')
%set(gca,'clim', [-3,7.5],'xlim',[-2,5]) %'clim', [-3,7.5]
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Primary Auditory')
subplot(3,1,3)
contourf(T,freq,tfV3(:,:,3),40,'linecolor','none')
%set(gca,'clim', [-3,7.5],'xlim',[-2,5]) %
colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Prefrontal Cortex')
suptitle('40Hz')


%% Saving data array
ID=file(12:16);
path1='D:\Results\2021_05_Pcdh9\Visual\SSVEP\';
IDs=strcat(path1,ID,'_SSVEP1.mat')
thisFile=cat(4,tfV1,tfV2,tfV3);
save(IDs, 'thisFile');

IDs1=strcat(path1,ID,'_SS1_ERP.mat')
thisFile1=ERP;
save(IDs1, 'thisFile1');
