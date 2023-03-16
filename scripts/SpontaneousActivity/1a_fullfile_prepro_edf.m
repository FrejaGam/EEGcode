%% Preprocessing script for .edf+
% preparing for scoring

[file,path]=uigetfile;

[header, signal_header, signal_cells] = blockEdfLoad(fullfile(path,file));

[ event_times, event_values, event_binary ] = parse_taini_edf_annotations( signal_header, signal_cells );

fs=signal_header(1).samples_in_record / header.data_record_duration;

%% Channels
% ch 13=PFC
EEG13 = cell2mat(signal_cells(1,13))';
EEG3 = cell2mat(signal_cells(1,3))'; 
EEG4 = cell2mat(signal_cells(1,4))';
EEG12 = cell2mat(signal_cells(1,12))';
EEG14 = cell2mat(signal_cells(1,14))';
% ch 8 & 11 both EMG
EMG8 = cell2mat(signal_cells(1,8))';
EMG11 = cell2mat(signal_cells(1,11))';

% all channels in one matrix
EEG=[EEG3;EEG4;EEG12;EEG13;EEG14;EMG8;EMG11];

% view of raw signal
figure(1)
hold on
subplot(2,1,1)
plot(EEG(4,:))
subplot(2,1,2)
plot(EEG(7,:))
hold off

%% filters and artifact rejection
EEGr=zeros(size(EEG));

for iCh = 1:size(EEG,1)

    %filter EEG data
    order = 4;
    hpFilt  = 0.5;
    lpFilt = 70;
    
    %Get rid of packeloss artificial low values and interpolate
    EEGCh = EEG(iCh,:);
    EEGCh (EEGCh < 6) = nan;                %ATTENTION thresholds may vary
    EEGCh (EEGCh > 10) = nan;
    Packloss(1,iCh) = (sum(isnan(EEGCh))/size(EEGCh,2)) *100;
    EEGres = resample(EEGCh,1:length(EEGCh));
    
    %Filter EEG data
    EEGHp = filterEEG(EEGres, order, hpFilt, fs, 'high', 2);
    EEGHLp = filterEEG(EEGHp, order, lpFilt, fs, 'low', 2);
    
    %reject artifacts
    average = mean(EEGHLp);
    stdev = std(EEGHLp);
    EEGHLp(EEGHLp < average - 3*stdev) = nan;
    EEGHLp(EEGHLp > average + 3*stdev) = nan;
    Packloss(2,iCh) = (sum(isnan(EEGHLp))/size(EEGHLp,2)) *100;
    EEGHLp = resample(EEGHLp,1:length(EEGHLp));
    
    EEGr(iCh,:)=EEGHLp;
end

%% plot artifact rejected signal

figure(2);
hold on
plot(EEG(1,:),'k')
plot(EEGr(1,:),'r')

%% Accuplot
% preview into the view of the signal

accuplot_small(EEGr,4,7,fs); % EEG matrix, ch EEG, ch EMG, fs

%% saving filtered signal
ID=file(12:16); % this line extracts the ID of the animal from the file name
path1='D:\Results\2020_11_Nrxn1\Visual\VEP\'; %define the path to save
IDs=strcat(path1,ID,'_VEP1.mat') %path+filename for the finished file
save(IDs, 'EEGr'); %this actually saves the data
