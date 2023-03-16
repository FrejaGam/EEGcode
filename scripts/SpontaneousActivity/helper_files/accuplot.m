function [G,fAxis,showFreqs]=accuplot(EEG, ch1,ch2, Mov, fs)
% function for creating the Accusleep function
% calls createSpectrogram.m, standardizeSR.m and processEMG.m (from Accusleep source code)
EEGr=EEG;
%% Accusleep
epochLen=1;
G = struct; % holds everything

G.originalSR = fs; % EEG/EMG sampling rate
G.SR = 128; % sampling rate used when calculating spectrogram and processed EMG
G.epochLen  = epochLen; % length of one epoch (spectrogram column) in seconds

G.EEG=EEGr(ch1,:); %EEGr(4,:) %state1(4,:); OPLET! which signal in which part of the matrix
G.EMG=EEGr(ch2,:); %EEGr(7,:)
G.mov=Mov;

% change variable names?
if length(G.EEG) ~= length(G.EMG)
    message = 'ERROR: EEG and EMG are different lengths';
    return
end

% create spectrogram and process EMG at a standard SR (128)
[spec, tAxis, fAxis] = createSpectrogram(standardizeSR(G.EEG, G.originalSR, G.SR), G.SR, G.epochLen);
G.processedEMG = processEMG(standardizeSR(G.EMG, G.originalSR, G.SR), G.SR, G.epochLen);
% set ceiling for EMG trace at 2.5 SD when plotting
G.cappedEMG = G.processedEMG;
emgCap = mean(G.cappedEMG) + 2.5*std(G.cappedEMG);
G.cappedEMG(G.cappedEMG > emgCap) = emgCap;

G.nbins = length(tAxis); % total number of time bins in the recording

% get spectrogram and time axes
showFreqs = find(fAxis <= 30); % only show frequencies under 30 Hz
G.specTs = (1:G.nbins)*G.epochLen - G.epochLen/2; % spectrogram time axis, in seconds
G.specTh = G.specTs./3600; % spectrogram time axis, in hours
G.spectrogram = spec(:,showFreqs); % our EEG spectrogram

%% Scaling
% take a sample of the spectrogram to help initialize the colormap
sampleBins = randperm(G.nbins, round(G.nbins/10));
specSample = reshape(spec(sampleBins,showFreqs),1,length(sampleBins)*length(showFreqs));
G.caxis1 = prctile(specSample,[2 98]); % percentile of signal used to scale the colorbar
G.cmax = G.caxis1(2);

% x scale for EEG that is not downsampled
dx=length(G.specTs)/length(G.EEG);
N=length(G.EEG);
G.tx=(0:N-1)*dx;

%% Figure
figure(1);
a1=subplot(411);
imagesc(G.specTs, fAxis(showFreqs), G.spectrogram',G.caxis1);
axis xy
a2=subplot(412);
hold on
plot(G.specTs,G.cappedEMG,'k');
set(gca,'YLim',[mean(G.cappedEMG)-1,mean(G.cappedEMG)+1],...
    'XLim',[0,max(G.specTs)]);
a3=subplot(413);
hold on
plot(Mov);       %in case of level
set(gca,'XLim',[0,max(length(Mov))]);
a4=subplot(414);
hold on
plot(G.tx,G.EEG);
set(gca,'XLim',[0,max(length(G.specTs))],'YLim',[-0.3,0.3]);

linkaxes([a1 a2 a3 a4],'x');
end