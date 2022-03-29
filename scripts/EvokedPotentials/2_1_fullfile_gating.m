%2020-02-27 saving section added
%2020-03-02 adapting to auditory
%2020-07-08 artifact rejection added

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

%% Draw triggers/ event points
%Convert event_times from seconds to datapoint by multiplying with the
%sampling frequency. Select every second points to select the 3's.
%transpose the array and remove first and last trigger. Round to end up at
%specific measurement points.
fs = signal_header(1).samples_in_record / header.data_record_duration;

event_times([1]) = []; %first trigger is the beginning of the recording
event_points = event_times * fs; 
event_points = event_points(:, 1:end)'; %every odd index
event_points = round(event_points);
plot(event_points, 1, 'r*')
hold on
ylim([0.5 1.5])
plot(event_points, 0.99, 'b*')
for ii = 1:length(event_points)
    text(event_points(ii),0.97,num2str(ii),'Color','k')
end

%% Divide over different ISIs
trig1 = event_points(1:2:199,:);
trig2 = event_points(201:2:399,:);
trig3 = event_points(401:2:599,:);
trig4 = event_points(601:2:799,:);
trig5 = event_points(801:2:1000,:);
trig6 = event_points(1001:2:1200,:);
trig7 = event_points(1201:2:1400,:);

%% Combining channels into one matrix
%select the EEG right channels from the struct and downsample to rSfreq
%3=Vright, 4=auditory, 12=auditory, 13=PFC , 14=Vleft (8 & 11= Muscle)
Ch3 = cell2mat(signal_cells(1,3))'; 
Ch4 = cell2mat(signal_cells(1,4))';
Ch12 = cell2mat(signal_cells(1,12))';
Ch13 = cell2mat(signal_cells(1,13))';
Ch14 = cell2mat(signal_cells(1,14))';
Ch8 = cell2mat(signal_cells(1,8))';
Ch11 = cell2mat(signal_cells(1,11))';

EEG = [Ch3;Ch4;Ch12;Ch13;Ch14;Ch8;Ch11];

figure(2)
hold on
subplot(2,2,1);
plot(EEG(2,:));
subplot(2,2,2);
plot(EEG(5,:));
subplot(2,2,3);
plot(EEG(1,:));
subplot(2,2,4);
plot(EEG(6,:));

%% Cut around triggers
% current standard i from -80 to 240
stimMat1 = bsxfun(@plus,trig1,-100:1500); 
stimMat2 = bsxfun(@plus,trig2,-100:1500);
stimMat3 = bsxfun(@plus,trig3,-100:1500); 
stimMat4 = bsxfun(@plus,trig4,-100:1500);
stimMat5 = bsxfun(@plus,trig5,-100:1500);
stimMat6 = bsxfun(@plus,trig6,-100:1500);
stimMat7 = bsxfun(@plus,trig7,-100:1500);

stimMat = cat(3,stimMat1, stimMat2,stimMat3,stimMat4,stimMat5,stimMat6,stimMat7);

T = (-100*1/rSfreq):1/rSfreq:(1500*1/rSfreq);

ISI1 = [1,2,3,4,5,7.5,10];

%% Filtering and artifact rejection
for iCh = 1:size(EEG,1)

    for iVolt = 1: size(stimMat,3)
        S2 = ISI1*108; %startpoint stimulus2
    
        %filter EEG data
        order = 4;
        hpFilt  = 0.5;
        lpFilt = 70;    %MEGs look awful if lowpass filter is increased
    
    
        %Get rid of packeloss artificial low values and interpolate
        EEGCh = EEG(iCh,:);
        EEGCh (EEGCh < 6) = nan;
        EEGCh (EEGCh > 10) = nan;
        tmp = EEGCh(stimMat(:,:,iVolt));  % make temp ERP to get trials which contain too long streches of package loss 
  
        for P = 1:size(tmp,1)
            temp = regionprops(isnan(tmp(P,[80:235,80+S2:235+S2])), 'area'); % find consecutive areas
            %Store max consecutive area per trial
            cons = max([(temp.Area)]);
            if cons > 1
            artP(P) = cons;
            else
            artP(P) = 0;
            end
        end
    
    % Find trials that have too long streches
    ArtP = find(artP>8);
    RemvP(iVolt,iCh) = length(ArtP);
    
    %Resample EEG data
    EEGres = resample(EEGCh,1:length(EEGCh));
    
    %Filter EEG data
    EEGHp = filterEEG(EEGres, order, hpFilt, rSfreq, 'high', 2);
    EEGHLp = filterEEG(EEGHp, order, lpFilt, rSfreq, 'low', 2);
    
            
    %reject artifacts
    tmp = EEGHLp(stimMat(:,:,iVolt));
    
    up = mean(tmp) + 5*std(tmp);
    down = mean(tmp) - 5*std(tmp);
    artA = tmp(:,[80:235,80+S2:235+S2])> up(:,[80:235,80+S2:235+S2]) | tmp(:,[80:235,80+S2:235+S2])<down(:,[80:235,80+S2:235+S2]);
    
    % Find trials that have artefacts
    [ArtA,col] = find(artA);
    RemvA(iVolt,iCh) = length(unique(ArtA));
    
    %Combine all to be removed trials
    Artt = unique([ArtA;ArtP']);
    ArtT(iVolt, iCh) = length(Artt);
    
    %Cut around triggers, average and baseline correct
    
        erp = EEGHLp(stimMat(:,:,iVolt));
        erp(Artt,:) = []; %Artt
        erp = mean(erp);
        erp = erp - mean(erp(1:78));
        ERP(iVolt,:,iCh) = erp;
        
        clear artA ArtA artP artP Artt
    end
    
end

%% visualization of raw EEG and artifact rejected EEG
figure(3)
hold on
plot(EEGres,'k')
plot(EEGHLp,'r')

%% plot triggers as subplots
figure(4)
hold on
subplot(7,1,1);
plot(T, ERP(1,:,3), 'k');
subplot(7,1,2);
plot(T, ERP(2,:,3), 'k');
subplot(7,1,3);
plot(T, ERP(3,:,3), 'k');
subplot(7,1,4);
plot(T, ERP(4,:,3), 'k');
subplot(7,1,5);
plot(T, ERP(5,:,3), 'k');
subplot(7,1,6);
plot(T, ERP(6,:,3), 'k');
subplot(7,1,7);
plot(T, ERP(7,:,3), 'k');

%% saving data matrix

ID=file(12:16);
path1='D:\Results\2021_05_Pcdh9\Auditory\gating\';
IDs=strcat(path1,ID,'_gating2.mat');
thisFile=ERP;
save(IDs, 'thisFile');

ArtT
