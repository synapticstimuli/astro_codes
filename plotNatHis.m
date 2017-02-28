function [natGrey, natStim] = plotNatHis(data,NatImages,Gratings)

% This funcion will load the NatImages and Gratings mat files, extract the
% blank screen responses and compare them to stimulus response for each
% protocol using a histogram plot. -RG 2016

%% Load mat files
clc
if nargin==0
    fprintf('Load data.mat file...\n')
    [fname, pathname] = uigetfile('.mat','Load data.mat file...');
    cd(pathname);
    
    try
        fprintf('Loading files...\n');
        natFile = strrep(fname,'_data','_NatImages');
        gratFile = strrep(fname,'_data','_Gratings');
        load (fname);
        load (natFile);
        load (gratFile);
        %         displayData(data);
    catch
        fprintf('Files not found, select another folder...\n');
    end
elseif nargin==3
    if ~exist('data','var')
        fprintf('Data var not loaded...\n');
        natFile = uigetfile('.mat','Loat data.mat...\n');
    else
        fname = strrep(data.filename,'.tif','_data.tif');
    end
    if ~exist('NatImages','var')
        fprintf('NatImages var not loaded...\n');
        natFile = uigetfile('.mat','Loat NatImages.mat...\n');
    else
        fprintf('NatImages var loaded...\n');
    end
    if ~exist('Gratings','var')
        fprintf('Gratings var not loaded...\n');
        natFile = uigetfile('.mat','Loat Gratings...\n');
    else
        fprintf('Gratings var loaded...\n');
    end
end

%% Read protocol from file name

% protInx = regexp(fname,'NS-');
% if strcmp(fname(protInx:protInx+2),'NS-')
%     data.stimDur    =   984;
%     nsProtocol        =   'NS1';
% elseif strcmp(fname(protInx:protInx+2),'NS2')
%     data.stimDur    =   984;
%     nsProtocol        =   'NS2';
% else
%     fprintf('Cannot determine protocol...\n');
%     nsProtocol = input('Enter protocol abbv (NS or NS2): \n','s');
% end

%% Input params
nsProtocol = 'NS1';
saveFlag        =   1;
plotFlag        =   0;
natFlag         =   0;

% Stim params
switch nsProtocol
    case 'NS1'
        Params.durBlank        = 4;
        Params.durBlank_Grats  = 6;
        Params.durStim         = 13;
        Params.durGratings     = 2;
        Params.repGrat         = 6;     % # of gratings stim repeats (total G's in V of NS stim protocol.
        %fprintf('Case = Astro\t\tTotal Time(s) = %d\n', 6.*((4*17)+(8*12)) );
    case 'NS2'
        Params.durBlank        = 4;
        Params.durBlank_Grats  = 1;
        Params.durStim         = 13;
        Params.durGratings     = 1;
        Params.repGrat         = 6;     % # of gratings stim repeats (total G's in V of NS stim protocol.
        %fprintf('Case = Astro\t\tTotal Time(s) = %d\n', 6.*((4*17)+(8*12)) );
end

fprintf('Protocol: %s\n', nsProtocol);
pause(3);
%% HouseKeeping
TotalDurStim = data.stimDur;

NumSamplsNew = TotalDurStim.*data.newFrameRate;

TimeEachTrial   =  (4*17)+(8*12);
SamplsEachTrial = TimeEachTrial*data.newFrameRate;
NumTrials       = NumSamplsNew./SamplsEachTrial;

%SamplsMovie     = 4*17*data.newFrameRate;
%SamplsGratings  = 8*12*data.newFrameRate;

NumRepsMov  = 4;

%% Partition nat movie responses
if natFlag
    
onDur = Params.durStim * data.newFrameRate;
offDur = Params.durBlank * data.newFrameRate; 
[dim1, dim2, dim3] = size(NatImages{1});

tic
ct = 0;

natGrey = [];
natStim = [];
for n = 1:NumTrials
    for i = 1:dim1
        for j = 1:dim2
            ct = ct+1;
            clc
            fprintf('Partitioning NS data...(%d of %d)\n', n, NumTrials);
            fprintf('Completed:%0.3f%% (%0.2fs)\n',(ct/((dim2^2)*NumTrials))*100,toc);
            natImgx = squeeze(NatImages{n}(i,j,:));
            repeats = reshape (natImgx, [(onDur + offDur) , NumRepsMov]);
            for tr=1:NumRepsMov
                foo = squeeze(repeats(:,tr));
                natG{tr}(i,j,:) = foo(onDur+1:onDur+offDur);
                natS{tr}(i,j,:) = foo(1:onDur);
                
            end
            
        end
        
    end
    
    natGx{n} = [];
    natSx{n} = [];
    for tr=1:NumRepsMov
        natGx{n} = cat(3, natGx{n},natG{tr});
        natSx{n} = cat(3, natSx{n},natS{tr});
    end

natStim = cat(3,natStim, natSx{n});
natGrey = cat(3,natGrey, natGx{n});
end

fprintf('Finished partitioning NS data...\n');
pause(3);
clearvars onDur offDur dim1 dim2 dim3 repeats

if saveFlag
    fprintf('Saving...\n')
    save (natFile, 'natStim', 'natGrey', '-append');
    fprintf('Done.\n');
end
end
%% Partition grating  responses
onDur = Params.durGratings * data.newFrameRate;
offDur = Params.durBlank_Grats * data.newFrameRate; 
[dim1, dim2, dim3] = size(Gratings{1});

tic
ct = 0;

gratGrey = [];
gratStim = [];
for n = 1:NumTrials
    for i = 1:dim1
        for j = 1:dim2
            ct = ct+1;
            clc
            fprintf('Partitioning gratings data...(%d of %d)\n', n, NumTrials);
            fprintf('Completed:%0.3f%% (%0.2fs)\n',(ct/((dim2^2)*NumTrials))*100,toc);
            gratImgx = squeeze(Gratings{n}(i,j,:));
            repeats = reshape (gratImgx, [(onDur + offDur) , NumRepsMov]);
            for tr=1:NumRepsMov
                foo = squeeze(repeats(:,tr));
                gratG{tr}(i,j,:) = foo(onDur+1:onDur+offDur);
                gratS{tr}(i,j,:) = foo(1:onDur);
            end
            
        end
        
    end
    
    gratGx{n} = [];
    gratSx{n} = [];
    for tr=1:NumRepsMov
        gratGx{n} = cat(3, gratGx{n},gratG{tr});
        gratSx{n} = cat(3, gratSx{n},gratS{tr});
    end

gratStim = cat(3,gratStim, gratSx{n});
gratGrey = cat(3,gratGrey, gratGx{n});
end

fprintf('Finished partitioning gratings data...\n');
pause(3);
clearvars onDur offDur dim1 dim2 dim3 repeats

if saveFlag
    fprintf('Saving...\n');
    save (natFile, 'gratStim', 'gratGrey', '-append');
    fprintf('Done.\n');
end


%% Plot max projection of on and off responses

if plotFlag
    %....Nat Stim On
    medNatStim=median(median(median(natStim)));
    natStim(natStim<medNatStim) = NaN;
    minvNatStim = min(min(min(natStim)));
    maxvNatStim = max(max(max(natStim)));
    natStim(isnan(natStim)) = minvNatStim-((maxvNatStim-minvNatStim)/5);
    
    medNatGrey=median(median(median(natGrey)));
    natGrey(natGrey<medNatGrey) = NaN;
    minvNatGrey = min(min(min(natGrey)));
    maxvNatGrey = max(max(max(natGrey)));
    natGrey(isnan(natGrey)) = minvNatGrey-((maxvNatGrey-minvNatGrey)/5);
    
    %  %RelGrat
    %  meanRelGrat=mean(mean(RelGrat));
    %  RelGrat(RelGrat<meanRelGrat) = NaN;
    %  minvGrat = min(min(RelGrat));
    %  maxvGrat = max(max(RelGrat));
    %  RelGrat(isnan(RelGrat)) = minvGrat-((maxvGrat-minvGrat)/5);
    
    h1=figure; set(h1,'color','w','position',[100 100 1000 400]);
    subplot(2,2,1);
    imagesc( vizMaxProjection(natStim) ); axis square; axis off; title('Stim on');
    subplot(2,2,2);
    imagesc( vizMaxProjection(natGrey) ); axis square; axis off; title('Stim Off');
    % subplot(2,2,3);
    % imagesc( vizMaxProjection(gratStim) ); axis square; axis off; title('Stim on');
    % subplot(2,2,4);
    % imagesc( vizMaxProjection(gratGrey) ); axis square; axis off; title('Stim Off');
end
