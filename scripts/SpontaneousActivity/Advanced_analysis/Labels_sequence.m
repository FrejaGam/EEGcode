%% Sequence Labels
% 2024-07-18
% this script asks for folder of inputs, and counts the different pairs of
% labels - the use is to quantify transitions as well as consecutive
% periods of the various states.

baseDir=uigetdir(pwd);
mat=dir(fullfile(baseDir,'*states.mat'));
nfiles=(length(mat));

sequences=cell(1,nfiles);

% load data from the different files and join the lables into pairs
for iFile=1:nfiles
    
    data=load(fullfile(baseDir,mat(iFile).name));
    
    labels=data.S.labels;
    
    Co=cat(2,labels(1:end-1,2),labels(2:end,2));
    Cojoin=sscanf(sprintf('%d%d,',[Co(:,1).';Co(:,2).']),'%d,');
    
    [GC,GR]=groupcounts(Cojoin);
    
    pairs=cat(2,GR,GC);
    
    sequences(iFile)={pairs};
    
end   

%% Combines a list of all possible pairs with the frequency of each
% all combinations of labels 8*8 [N a1 a2 a3 REM nREM1 nREM2 nREM3]
v=[78 97 98 99 108 109 110 111]; 

c1=nchoosek(v,2);
c2=cat(2,c1(:,2),c1(:,1));
combo=cat(1,c1,c2,cat(1,v,v)'); 
stand=sscanf(sprintf('%d%d,',[combo(:,1).';combo(:,2).']),'%d,');
stand=sort(stand);

totSeq=[];

for i=1:size(sequences,2)
    x=cell2mat(sequences(1,i));
    
    found=ismember(stand,x(:,1));
    result(found)=x(:,2);
    GC=result';
    
    totSeq=cat(2,totSeq,GC);
end

figure(1)
bar(mean(totSeq,2))

%% Saving the file
path1='/Users/kaslab/Desktop/Test files Pcdh9/';
IDs=strcat(path1,'Pcdh9_seqLabels.mat');
save(IDs,'totSeq');