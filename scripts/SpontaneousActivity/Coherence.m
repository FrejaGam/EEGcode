%% Coherence from scored data
% 2024-07-18
% loops over a folder of data that has been distributed into vigilance
% states

% requires no extra functions

baseDir=uigetdir(pwd); % asks for input folder
mat=dir(fullfile(baseDir,'*.mat')); % finds .mat files in that folder
nfiles=length(mat); % number of relevant files

fs=1084.7; %sampling frequency in Hz

for thisFile=1:nfiles
    in_mat=fullfile(baseDir,mat(thisFile).name);
    data=load(in_mat);
%%
    a1=cell2mat(data.S.states(1,1));
    a2=cell2mat(data.S.states(1,2));
    a3=cell2mat(data.S.states(1,3));
    rem=cell2mat(data.S.states(1,4));
    nrem1=cell2mat(data.S.states(1,5));
    nrem2=cell2mat(data.S.states(1,6));
    nrem3=cell2mat(data.S.states(1,7));



%% Coherence looping through all states
    l=200000; %set by computing limitations
    Coh=[];

    list={a1,a2,a3,rem,nrem1,nrem2,nrem3}; %relevant states
    for i=1:length(list)
        raw=list{i};
        Coh=[];
        if size(raw,2)<l
            l=size(raw,2);
        else
            l=200000;
        end
        % relevant pairs should be built by combining i1 and i2
        for i1=[1 2]
            for i2=[5 3 2]
                x=raw(i1,1:l);
                y=raw(i2,1:l);
                [wcoh,wcs,f]=wcoherence(x,y,fs);
                Coh_m=mean(wcoh,2);
                Coh=cat(2,Coh,Coh_m);
            end
        end
        cohs(i)={Coh};
    end 

    %% Standardizing outcome coherence
    %cohs is a struct of fields, that might be inequal in length, depending
    %on the amount of data in each state
    M=min(cellfun('size',cohs,1));
    TA=[];
    for i=1:length(list)% [1,2,4,5,6,7]
        raw=cell2mat(cohs(1,i));
        TA=cat(2,TA,raw(1:M,:)); %depends on cohs
    end

    Coh=TA;
    %% insert Nan ONLY if a state is missing
    
    %repl = nan(177,6);
    %Coh=cat(2,Coh(:,1:12),repl,Coh(:,13:end));
    
    %% Saving
    ID=in_mat(53:64); % relevant part of specific file name
    path1='/Users/kaslab/Desktop/Test files Pcdh9/Pcdh9_coherence/'; % path to output folder
    IDs=strcat(path1,ID,'_coh.mat'); % name of new file
    save(IDs, 'Coh'); % actual saving function
end