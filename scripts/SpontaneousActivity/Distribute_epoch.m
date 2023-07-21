%% assigning scored data from labels
%% Load
close all
clear all

baseDir=uigetdir(pwd);      % find labels
mat=dir(fullfile(baseDir,'*'));
nfiles=(length(mat));

for iFile=3:length(mat)
    data=load(fullfile(baseDir,mat(iFile).name));
    m=data.S.labels;
    IS=data.S.sleep;

    id=mat(iFile).name(3:8);
    fn=strcat('*',id,'*');
    Gdir='D:\Sleep_scoring\Nrxn1_23\singleFiles\';
    spec=dir(fullfile(Gdir,fn));
    EEG=load(fullfile(Gdir,spec.name));
    EEGr=EEG.aE.EEG;
    fs=454; %
    %% data for each state
    labels=[97:99,108:111]; %108:111
    states=cell(1,8);

    for i=1:7 %labels of states
        indV=find(m(:,2)==labels(i));
        state=[];
        for ij=1:size(indV) %each epoch with a specific label
            ee=round(indV(ij)*fs); %epoch end
            es=round((indV(ij)-1)*fs); %epoch start
            if indV(ij)==1
                epoch=EEGr(:,1:ee);
            else
                epoch=EEGr(:,es:ee);
            end
        state=cat(2,state,epoch);
        end
        states(i)={state};
    end

    for ij=1:size(IS,2) %each epoch with a specific label
        ee=round(IS(ij)*fs); %epoch end
        es=round((IS(ij)-1)*fs); %epoch start
        if IS(ij)==1
           epoch=EEGr(:,1:ee);
        else
           epoch=EEGr(:,es:ee);
        end
        state=cat(2,state,epoch);
    end
    states(8)={state};
    %% saving
    S.states=states;
    S.labels=m;

    path1='D:\Sleep_scoring\Nrxn1_23\scored\';
    IDs=strcat(path1,'2_',id,'_states.mat');
    save(IDs, 'S');
end