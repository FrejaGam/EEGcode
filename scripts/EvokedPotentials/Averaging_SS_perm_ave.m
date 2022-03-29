%% load color maps
clear all

baseDir=uigetdir(pwd);
cd(baseDir);

mat = dir('*SSVEP3_col*');
 for q = 1:length(mat) 
     data(q)=load(mat(q).name);     
 end
 
min_freq =  5; % in Hz
max_freq = 70; % in HZ
num_freq = 65; % in count
freq = linspace(min_freq,max_freq,num_freq);
 
rSfreq = 1085;
T = (-3500*1/rSfreq):1/rSfreq:(4500*1/rSfreq);
 
cd('D:\MATLAB')
 
%% Figure
i=1;
 
figure(1)
subplot(311)
contourf(T,freq,data(1).S.rHET(:,:,i),40,'linecolor','none');
subplot(312)
contourf(T,freq,data(2).S.rHET(:,:,i),40,'linecolor','none');
subplot(313)
contourf(T,freq,data(3).S.rHET(:,:,i),40,'linecolor','none');

%% mean HET
HET=cat(4,data(1).S.rHET,data(2).S.rHET,data(3).S.rHET);
mHET=mean(HET,4);

figure(2)
subplot(311)
contourf(T,freq,mHET(:,:,1),40,'linecolor','none');
subplot(312)
contourf(T,freq,mHET(:,:,2),40,'linecolor','none');
subplot(313)
contourf(T,freq,mHET(:,:,3),40,'linecolor','none');

%% mean KO
KO=cat(4,data(1).S.rKO,data(2).S.rKO,data(3).S.rKO);
mKO=mean(KO,4);

figure(2)
subplot(311)
contourf(T,freq,mKO(:,:,1),40,'linecolor','none');
subplot(312)
contourf(T,freq,mKO(:,:,2),40,'linecolor','none');
subplot(313)
contourf(T,freq,mKO(:,:,3),40,'linecolor','none');

%% cluster-based permutation stat HET
fHET=cat(4,data(1).S.mFakeHET,data(2).S.mFakeHET,data(3).S.mFakeHET);
fHETm=squeeze(mean(fHET,4));
fHETs=cat(4,data(1).S.sdFakeHET,data(2).S.sdFakeHET,data(3).S.sdFakeHET);
fHETsd=squeeze(mean(fHETs,4));

sigThresh=norminv(1-0.05/2);
% zmap
[zmap,zthreshMap]=deal(mHET-fHETm ./ fHETsd);
zthreshMap(abs(zmap)<sigThresh)=0;

% i is condition
i=3;
islands=bwconncomp(logical(zthreshMap(:,:,i)));

figure(3)
hold on
contourf(T,freq,zmap(:,:,i),40,'linecolor','none')
contour(T,freq,zthreshMap(:,:,i),1,'linecolor','k')
caxis([-3 3])
%colorbar

%% cluster-based permutation stat KO
fKO=cat(4,data(1).S.mFakeKO,data(2).S.mFakeKO,data(3).S.mFakeKO);
fKOm=squeeze(mean(fKO,4));
fKOs=cat(4,data(1).S.sdFakeKO,data(2).S.sdFakeKO,data(3).S.sdFakeKO);
fKOsd=squeeze(mean(fKOs,4));

% zmap
[zmap_ko,zthreshMap_ko]=deal(mKO-fKOm ./ fKOsd);
zthreshMap_ko(abs(zmap_ko)<sigThresh)=0;

% i is condition
ii=3;
islands_ko=bwconncomp(logical(zthreshMap_ko(:,:,ii)));

figure(4)
hold on
contourf(T,freq,zmap_ko(:,:,ii),40,'linecolor','none')
contour(T,freq,zthreshMap_ko(:,:,ii),1,'linecolor','k')
caxis([-3 3])
%colorbar