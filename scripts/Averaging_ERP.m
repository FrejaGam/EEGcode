%% Averaging SingleERP
clear all
close all

baseDir=uigetdir(pwd);
cd(baseDir);

fs = 1085;
T = (-80*1/fs):1/fs:(240*1/fs);
 
 mat = dir('*VEP4r.mat');
 for q = 1:length(mat) 
     data(q)=load(mat(q).name);     
 end

cd('D:\MATLAB')  
%% Divide into groups and make grand averages
 %3=Vleft, 4=Aleft, 12=Aright, 13=PFC , 14=Vright (8 & 11= Muscle)
 WT=cat(3,data(7).thisFile(:,1),data(7).thisFile(:,5),data(13).thisFile(:,1),data(13).thisFile(:,5), ...
    data(12).thisFile(:,1),data(12).thisFile(:,5),data(14).thisFile(:,1),data(14).thisFile(:,5), ...
    data(10).thisFile(:,1),data(10).thisFile(:,5),data(1).thisFile(:,1),data(1).thisFile(:,5));
meanWT=mean(WT,3);
sdWT=std(WT,0,3);
 
 KO=cat(3,data(16).thisFile(:,1),data(16).thisFile(:,5),data(15).thisFile(:,1),data(15).thisFile(:,5), ...
     data(11).thisFile(:,1),data(11).thisFile(:,5),data(19).thisFile(:,1),data(19).thisFile(:,5), ...
     data(17).thisFile(:,1),data(17).thisFile(:,5),data(18).thisFile(:,1),data(18).thisFile(:,5));
meanKO=mean(KO,3);
sdKO=std(KO,0,3);

HET=cat(3,data(2).thisFile(:,1),data(2).thisFile(:,5),data(8).thisFile(:,1),data(8).thisFile(:,5), ...
    data(9).thisFile(:,1),data(9).thisFile(:,5),data(3).thisFile(:,1),data(3).thisFile(:,5), ...
    data(4).thisFile(:,1),data(4).thisFile(:,5),data(5).thisFile(:,1),data(5).thisFile(:,5), ...
    data(6).thisFile(:,1),data(6).thisFile(:,5));
meanHet=mean(HET,3);
sdHet=std(HET,0,3);

figure(1)
hold on
plot(T,meanWT,'b')
plot(T,meanKO,'k')
plot(T,meanHet,'r')
title('Grand averages: VEP','FontSize',20)
legend({'WT','Hom KO','Het'},'FontSize',15)

figure(2)
hold on
shadedErrorBar(T,meanWT,sdWT,'lineprops','-b')
shadedErrorBar(T,meanKO,sdKO,'lineprops','-k')
shadedErrorBar(T,meanHet,sdHet,'lineprops','-r')
title('Grand averages: VEP','FontSize',20)
legend({'WT','Hom KO','Het'},'FontSize',15)

%% Checking quality of channels and time windows

for j = 1:size(WT,3)
    figure(j)
    hold on
    plot(T,WT(:,:,j))
    xline(T(80))
    xline(T(90))
    xline(T(97))
    xline(T(120))
    xline(T(127))
    xline(T(140))
    xline(T(250))
    hold off
end

%% Peaks definitions
gr=cat(3,WT,KO,HET); 

peaksA = zeros(size(gr,3),7);
peaksI = zeros(size(gr,3),7);

for iAni = 1:size(gr,3)
     
     [P1a,P1i] = max(gr(80:90,:,iAni),[],1);
     [N1a,N1i] = min(gr(91:97,:,iAni),[],1);
     [P2a,P2i] = max(gr(98:120,:,iAni),[],1);
     [N2a,N2i] = min(gr(98:120,:,iAni),[],1);
     [P3a,P3i] = max(gr(121:127,:,iAni),[],1);
     [N3a,N3i] = min(gr(128:140,:,iAni),[],1);
     [P4a,P4i] = max(gr(141:250,:,iAni),[],1);
     
     peaksA(iAni,1)=P1a;
     peaksA(iAni,2)=N1a;
     peaksA(iAni,3)=P2a;
     peaksA(iAni,4)=N2a;
     peaksA(iAni,5)=P3a;
     peaksA(iAni,6)=N3a;
     peaksA(iAni,7)=P4a;
     
     peaksI(iAni,1)=P1i;
     peaksI(iAni,2)=N1i;
     peaksI(iAni,3)=P2i;
     peaksI(iAni,4)=N2i;
     peaksI(iAni,5)=P3i;
     peaksI(iAni,6)=N3i;
     peaksI(iAni,7)=P4i;
end     
