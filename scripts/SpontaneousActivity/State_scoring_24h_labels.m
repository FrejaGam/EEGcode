%% State scoring
% 23-01-10 updated with frequency band requirement for REM (FGO)
%% load section
close all
clear all

[file,path]=uigetfile;
in_mat=fullfile(path,file);
data=load(in_mat);

EEGr=data.aE.EEG; %data.aE.EEG data.EEGr
Mov=data.aE.Mov; %data.aE.Mov mov=0 (if no mov)
fs=454; %1028 1084.7

[G,fAxis,showFreqs]=accuplot(EEGr,4,7,Mov,fs); % EEG, ch EEG, ch EMG, mov sensor, fs
%% Mean of spectrogram data
% and crude artifact rejection
m=mean(G.spectrogram,2);

for i=1:size(G.spectrogram,1)
    if m(i,1)>0.002
        m(i,2)='N';
    elseif m(i,1)<0.000002
        m(i,2)='N';
    else
        m(i,2)='a';
    end
end

figure(2)
plot(m(:,1))
%% smooth EMG & Find transitions
trans=6; % estimate or count the number of transitions

emgCap = -6;
G.cappedEMG(G.cappedEMG < emgCap) = emgCap;

figure(3)
findchangepts(m(:,1).*G.cappedEMG','MaxNumChanges',trans)

figure(4)
findchangepts(G.cappedEMG,'MaxNumChanges',trans) %less affected by noise

% changepoints found in either EMG or EEG
cpEMG=findchangepts(G.cappedEMG','MaxNumChanges',trans);
cpEEG=findchangepts(m(:,1).*G.cappedEMG','MaxNumChanges',trans);


%% divide sections into active and sleep
% get rid of clusters of change point mean & round
% use either changepoints in cpEMG or cpEEG or use epoch numbers directly
IA=[1:716,2373:4199,4666:4889,6829:7075,8454:9229,10018:10504,...
    12517:13176,15370:15636,17468:18118,18647:19269,20108:21010,...
    24462:size(m,1)];    %cpEEG(1) size(m,1)
Act=G.spectrogram(IA,:);
IS=setdiff(1:size(G.spectrogram,1),IA);
Sle=G.spectrogram(IS,:);

figure(5)
subplot(211)
imagesc(1:size(Act,1),fAxis(showFreqs), Act',G.caxis1)
axis xy
subplot(212)
imagesc(1:size(Sle,1), fAxis(showFreqs), Sle',G.caxis1)
axis xy

%pay attention to sleep in the active state

%% Scoring active state
% a=97 awake/active b=98 awake c=99 ressting state

for i=1:size(Act,1)
    if m(IA(i),1)<0.00012 && m(IA(i),1)>0.000002 %check figure 2
        bin1=mean(Act(i,1:15)); %low delta, or delta 0.5-4 Hz [1:21]
        bin2=mean(Act(i,16:30)); %broad theta, or theta 4-8 Hz [21:41]
        bin3=mean(Act(i,31:60)); %alpha 8-12 Hz  [41:61]
        if (bin1>bin2) && (bin1>bin3)
            m(IA(i),2)='a';
        elseif bin2>bin1 && bin2>bin3
            m(IA(i),2)='b';
        elseif bin3>bin2 && bin3>bin1
            m(IA(i),2)='c';
        else
            m(IA(i),2)='U';
        end
    else
        m(IA(i),2)='N';
    end
end

% plot of ASCII labels from the loop above
figure(6)
plot(m(:,2))

%% sleep scoring
% When all the awake states have been scored it is time to look at sleep
ms=mean(Sle,2);

figure(7)
plot(ms)
%% sleep labelling
%REM sleep l=108, slow-wave m=109, intermediate n=110, fast o=111

for i=1:size(Sle,1)
    if ms(i)<0.00025 && ms(i)>0.000002 %check upper bound especially
        bin1=mean(Sle(i,1:60)); %1:60
        binC=mean(Sle(i,100:120)); %100:120
        binR=mean(Sle(i,30:50));
        binRC=mean(Sle(i,1:30));
        if ms(i)<0.00007 && binR>binRC %check figure 7
            m(IS(i),2)='l';
        elseif bin1/binC>20 
            m(IS(i),2)='m';
        elseif bin1/binC<20 && bin1/binC>10 % bin1<binC && binR/binC<0.25
            m(IS(i),2)='n';
        elseif bin1/binC<10
            m(IS(i),2)='o';
        else
            m(IS(i),2)='U';
        end
    else
        m(IS(i),2)='N';
    end
end

figure(8)
plot(m(:,2))
            
%% State testing
% Here you can check three epochs of interest, does the labelling look
% reasonable
a=14082;
b=15362;
c=17073;

figure(9)
subplot(311)
plot(fAxis(showFreqs),G.spectrogram(a,:))%
subplot(312)
plot([0:0.2:30],G.spectrogram(b,:))
subplot(313)
plot(G.spectrogram(c,:))

%%
S.sleep=squeeze(IS); % Sleep in it's natural order
S.labels=m;
%% saving the list of labels
% afterwards the filtered data can be distributed according to the label by
% using the script Distribute_epoch.m
ID=file(8:13);
path1='D:\Sleep_scoring\Nrxn1_23\scored\labels\';
IDs=strcat(path1,'2_',ID,'_labs.mat');
save(IDs, 'S');

%% In case of video
%Load video and relevant frames
[fileV,pathV]=uigetfile;
v=VideoReader(fullfile(pathV,fileV));
% load .csv with frame times from Spike2
[fileT,pathT]=uigetfile;
txt=csvread(fullfile(pathT,fileT));
%% figure of interesting frames
%timing from Accusleep (figure 1)
J=1147; %timepoint of interest
[c,indV]=min(abs(txt(:,2)-J)); %csv: two columns; ch: one column
% read first frame, and the 19 following
frame=read(v,[indV indV+19]);            % (object, framenumber(s))
for i=1:size(frame,4)
    figure(15)
    subplot(4,5,i)
    image(frame(100:end,200:450,:,i)) %can be cropped here 
end

%alternative (doesn't adjust timing!!)
[row,col]=find(m==98); % ASCII labels

for i=1:20
    j=row(i);
    frame=read(v,j);
    figure(10)
    subplot(4,5,i)
    image(frame(100:end,250:450,:,:))
end

%% state spectrogram

for i=1:20
    R=[J:J+19];
    figure(11)
    subplot(4,5,i)
    plot(fAxis(1:100),G.spectrogram(R(i),1:100))
end

%% Indices of states (optional)
indV=find(m(:,2)~=78);
state=[];
for i=1:size(indV)
    epoch=EEGr(:,(indV(i)-1)*fs:indV(i)*fs);
    state=cat(2,state,epoch);
end
figure(12)
plot(state(1,:))