%% Permutations for SS response
clear all
close all

baseDir=uigetdir(pwd);
cd(baseDir);

rSfreq = 1085;
T = (-3500*1/rSfreq):1/rSfreq:(4500*1/rSfreq);
 
mat = dir('*SSVEP4*.mat');

for q = 1:length(mat) 
    data(q)=load(mat(q).name);     
end
 
min_freq =  5; % in Hz
max_freq = 70; % in HZ
num_freq = 65; % in count
freq = linspace(min_freq,max_freq,num_freq);
 
 cd('D:\MATLAB')
 
%% load all relevant data into one big array
 idata=[];
 for i=1:length(mat)
     idata=cat(3,idata,data(i).thisFile(:,:,1,:),data(i).thisFile(:,:,5,:));
 end
 
%% permutation to create a reference for later z-stats
 rng(1);            %random number generator for the samples function used later
 WT=[];             %relevant groups
 KO=[];
 HET=[];
 
 for perm=1:100
     samples=randperm(size(idata,3));
     
     WT=cat(3,WT,mean(idata(:,3000:end,samples(1:10),:),3));
     KO=cat(3,KO,mean(idata(:,3000:end,samples(11:18),:),3));
     HET=cat(3,HET,mean(idata(:,3000:end,samples(19:end),:),3));
 end

%% Differences and means for randomly assigned data (Fake)
Fake_HET=HET-WT;
Fake_KO=KO-WT;

S.mFakeHET=squeeze(mean(Fake_HET,3));
S.sdFakeHET=squeeze(std(Fake_HET,[],3));
S.mFakeKO=squeeze(mean(Fake_KO,3));
S.sdFakeKO=squeeze(std(Fake_KO,[],3));

%% threshold for cluster sizes

mFake=(Fake_HET+Fake_KO)/2;
m1Fake=mean(mFake);
sdFake=std(mFake);

clustersizes=zeros(100,1);
sigThresh=norminv(1-0.05/2);

for permi=1:100
    % z diff
    zdiffFake=(mFake(:,:,permi,:)-m1Fake)./sdFake;
    % threshold
    zdiffFake(abs(zdiffFake)<sigThresh)=0;
    % find clusters
    islands=bwconncomp(logical(zdiffFake));
    % record cluster sizes
    clustNs=cellfun(@length,islands.PixelIdxList);
    clustersizes(permi)=max(clustNs);
end

S.clustthresh=prctile(clustersizes,95);

figure(1)
hold on
histogram(clustersizes)
plot([1 1]*S.clustthresh,get(gca,'ylim'),'r','linew',3)

%% Real groups
%b1: WT [1:2,13:14,19:22,25:26,37:38] KO [23:24,27:36] HET [3:12,15:18]
%b2: WT [9:10,15:16,19:20] KO [1:2,7:8,11:14] HET [3:6,17:18,21:22]
%b3: WT [3:4,7:10,13:14,21:22] KO [1:2,5:6,11:12,15:16] HET [17:20]

MrWT=squeeze(mean(idata(:,3000:end,[3:4,7:10,13:14,21:22],:),3));
MrKO=squeeze(mean(idata(:,3000:end,[1:2,5:6,11:12,15:16],:),3));
MrHET=squeeze(mean(idata(:,3000:end,[17:20],:),3));

S.rKO=MrKO-MrWT;
S.rHET=MrHET-MrWT;


path='D:\Results\2020_04_Nrxn1\SSVEP_mean\new\';
IDs=strcat(path,'b3_SSVEP4_color_perm.mat');
save(IDs,'S');

%%
contourf(T(3000:end),freq,S.mFakeKO(:,:,3),40,'linecolor','none');
