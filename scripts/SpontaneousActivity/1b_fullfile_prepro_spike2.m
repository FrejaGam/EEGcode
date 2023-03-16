%% load of merged matfile from Spike2
% for resting state analysis with Accusleep

clear all
close all

[file,path]=uigetfile;
in_mat=fullfile(path,file);
data=load(in_mat);
f_names=fieldnames(data);

%% load each channel
%3=Vleft, 4=auditoryLeft, 12=audRight, 13=PFC , 14=Vright (8 & 11= Muscle)

[EEG_3, EEG_4, EEG_14, EEG_12, EEG_13]=deal([]);
[EMG_8, EMG_11, Mov]=deal([]);

EEG_3=data.(f_names{20});       
EEG_3=EEG_3.values;             
EEG_4=data.(f_names{21});
EEG_4=EEG_4.values;
EEG_14=data.(f_names{6});
EEG_14=EEG_14.values;
EEG_12=data.(f_names{4});
EEG_12=EEG_12.values;
EEG_13=data.(f_names{5});
EEG_13=EEG_13.values;
    
EMG_8=data.(f_names{26});
EMG_8=EMG_8.values;
EMG_11=data.(f_names{3});
EMG_11=EMG_11.values;


%% Movement sensor    
% function for data from motion sensors in the level format
Mov1=data.(f_names{13}); 

Mov=zeros(round(Mov1.times(end)),1);   
Mov(round(Mov1.times))=1;

%% EEG matrix 
fs=454;    %rounded sampling frequency
EEG = [EEG_3,EEG_4,EEG_12, EEG_13, EEG_14, EMG_8,EMG_11]; %matrix with all relevant channels
EEG=EEG';

%% Filter to Accusleep
EEGr=zeros(size(EEG));

for iCh = 1:size(EEG,1)

    %filter EEG data
    order = 4;
    hpFilt  = 0.5;
    lpFilt = 70;
    
    %Get rid of packeloss artificial low values and interpolate
    EEGCh = EEG(iCh,:);
    EEGCh (EEGCh < 7) = nan; 
    EEGCh (EEGCh > 8.8) = nan;                
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
plot(EEG(4,:),'k')
plot(EEGr(4,:),'r')

%% Accusleep
% preview of signal
[G,fAxis,showFreqs]=accuplot(EEGr,4,7,mov,fs); % EEG matrix, ch EEG, ch EMG, motion sensor, fs
%% Structure for saving all EEG, motion data and package loss
aE.EEG=EEGr;
aE.Mov=Mov;
aE.loss=Packloss;

%% Saving the EEG file
ID=file(1:12); %subject ID
path1='D:\Sleep_scoring\Nrxn1_23\singleFiles\'; % path for saving
IDs=strcat(path1,ID,'_restState.mat'); %concatenated name of new file
save(IDs, 'aE');
