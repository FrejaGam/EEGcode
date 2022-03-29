%% Averaging gating
clear all
close all

baseDir=uigetdir(pwd);
cd(baseDir);

rSfreq = 1085;
T = (-100*1/rSfreq):1/rSfreq:(1500*1/rSfreq);
 
 mat = dir('*gating1.mat');
 for q = 1:length(mat) 
     data(q)=load(mat(q).name);     
 end

cd('D:\MATLAB')  
%% Divide into groups and make grand averages
%3=Vleft, 4=Aleft, 12=Aright, 13=PFC , 14=Vright (8 & 11= Muscle)
%interstimulus interval (ISI) 1:5 for visual and 1:7 for auditory
i=5;
 
WT=cat(3,data(5).thisFile(i,:,1),data(5).thisFile(i,:,5),data(8).thisFile(i,:,1),data(8).thisFile(i,:,5), ...
    data(10).thisFile(i,:,1),data(10).thisFile(i,:,5));
meanWT=mean(WT,3);
sdWT=std(WT,0,3);
 
KO=cat(3,data(1).thisFile(i,:,1),data(1).thisFile(i,:,5),data(4).thisFile(i,:,1),data(4).thisFile(i,:,5), ...
     data(6).thisFile(i,:,1),data(6).thisFile(i,:,5),data(7).thisFile(i,:,1),data(7).thisFile(i,:,5));
meanKO=mean(KO,3);
sdKO=std(KO,0,3);

HET=cat(3,data(2).thisFile(i,:,1),data(2).thisFile(i,:,5),data(3).thisFile(i,:,1),data(3).thisFile(i,:,5), ...
    data(9).thisFile(i,:,1),data(9).thisFile(i,:,5),data(11).thisFile(i,:,1),data(11).thisFile(i,:,5));
meanHet=mean(HET,3);
sdHet=std(HET,0,3);

figure(1)
hold on
plot(T,meanWT,'b')
plot(T,meanKO,'k')
plot(T,meanHet,'r')
title('Grand averages: Gating','FontSize',20)
legend({'WT','Hom KO','Het'},'FontSize',15)

figure(2)
hold on
shadedErrorBar(T,meanWT,sdWT,'lineprops','-b')
shadedErrorBar(T,meanKO,sdKO,'lineprops','-k')
shadedErrorBar(T,meanHet,sdHet,'lineprops','-r')
title('Grand averages: Gating','FontSize',20)
legend({'WT','Hom KO','Het'},'FontSize',15)

%% Checking quality of channels

for j = 1:size(WT,3)
    figure(j)
    hold on
    plot(T,WT(:,:,j))
    hold off
end

%% Definition of responses
ISI = [0.1,0.2,0.3,0.4,0.5];
ISI = round(ISI(i) *rSfreq);
 
figure(11)
hold on
plot(T,meanWT,'b')
xline(T(80))
xline(T(180))
xline(T(ISI+160),'b')
xline(T(ISI+320),'b')

%% Calculation of gating
gr=cat(3,WT,KO,HET); 

SenGat = zeros(1,size(gr,3));
 
 for iAni = 1:size(gr,3)
     
     P1C = squeeze(max(gr(:,80:180,iAni),[],2));
     N1C = squeeze(min(gr(:,80:180,iAni),[],2));
     DiffC = P1C - N1C;
     P1N1C(iAni) = DiffC;
          
     P1T = squeeze(max(gr(:,ISI+160:ISI+320,iAni),[],2));
     N1T = squeeze(min(gr(:,ISI+160:ISI+320,iAni),[],2));
     DiffT = P1T - N1T;
     P1N1T(iAni) = DiffT;
       
     SenGat(iAni) = [P1N1T(iAni)./P1N1C(iAni)];
 end
 
