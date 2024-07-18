%% Averaging labels
% 2024-07-18
% This script averages labels from the scoring script

baseDir=uigetdir(pwd);
mat=dir(fullfile(baseDir,'*states.mat'));
nfiles=length(mat);

totGC=[];

for i=1:4:nfiles
    iGC=zeros(8,1);
    for j=0:3
        data=load(fullfile(baseDir,mat(i+j).name));
        labels=data.S.labels(:,2);
        [GC,GR]=groupcounts(labels);

        if length(GR)<8
            stand=[78;97;98;99;108;109;110;111];
            found=ismember(stand,GR);
            result(found)=GC;
            GC=result';
        end
        iGC=iGC+GC;
    end
    totGC=cat(2,totGC,iGC);
end


%% Normalization
norm_GC=totGC./sum(totGC,1);

%% Stacked plot of mean states
nm1=mean(norm_GC,2);

grMn=cat(2,nm1,nm1)'; %concatenate data from each group 
groups=categorical({'group 1' 'group 2'});

figure(2)
hold on
b=bar(groups',grMn,'stacked');
labels = { 'Noise', 'delta awake', 'theta awake',  'alpha awake', 'REM sleep', 'slow nREM', 'medium nREM', 'fast nREM' }; 
legend(labels,'location','eastoutside');
ypos = transpose(cat(1,b.YEndPoints)-[b(1).YEndPoints/2;diff(cat(1,b.YEndPoints))/2]);      % Calculate text positions
text(cat(2,b.XEndPoints),ypos(:),arrayfun(@(x) sprintf('%.2f',x),(cat(2,b.YData)),'uni',0),...
    'HorizontalAlignment','center','VerticalAlignment','middle');
