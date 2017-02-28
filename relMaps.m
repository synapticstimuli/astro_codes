function [data,NatImages,RelNat,Gratings] = relMaps(fname,protocol)
% Processes the data file and calculates the reliability of the responses
% to both NS and Gratings. Distuinguishes between different stim protocols 
% (eg, NS1 or NS2). Runs about 1.5 hrs per file. From RVR, modified for 
% astros by RG.

%% Load data file..................................................
clc
if nargin==0
    fprintf('Select data.mat file...\n')
    [fname, pathname] = uigetfile('*data.mat','Load .mat file containing data structure...');
    cd(pathname);
    fprintf('Loading file...\n');
    load (fname, 'data')
    protocol = [];
else
    if exist('data','var') == 1
        fprintf('Data structure loaded...');
        disp(data);
    elseif exist('data','var') == 0
        fprintf('Loading file...\n');
        load (fname, 'data')
        fprintf('Data structure loaded...\n');
%         disp(data);
    end
end
% User variables
displayFlag     =   0;
saveFlag        =   1;
plotFlag        =   1;
protocol        =   'NS';       % Default protocol.
%% Read protocol from file name....................................
% Reads file name and determines protocol - this will set stimulus parameters.

if ~isempty(strfind(fname,'-NS-'))
    data.stimDur    =   984;
    protocol        =   'NS1';
    fprintf('Protocol: %s\n', protocol);
    pause(3);
elseif ~isempty(strfind(fname,'-NS2-'))
    data.stimDur    =   984;
    protocol        =   'NS2';
    fprintf('Protocol: %s\n', protocol);
    pause(3);
elseif ~isempty(strfind(fname,'-NSK-'))
    data.stimDur    =   984;
    protocol        =   'NS1';
    fprintf('Protocol: %s\n', protocol);
    pause(3);
else
    fprintf('Cannot determine protocol...\n');
    fprintf('Running default protocol: %s...\n',protocol);
    fprintf('Protocol: %s\n', protocol);
    pause(3);
%     protocol = input('Enter protocol abbv (NS or NS2): \n','s');
end

%% Stimulus parameters.............................................
% Very important. Must match stimulus protocol parameters for scrip to
% funciton properly. Refer to NatScenes_Gratings(var) script(s).
switch protocol
    case 'NS1'
        Params.durBlank        = 4;
        Params.durBlank_Grats  = 6;
        Params.durStim         = 13;
        Params.durGratings     = 2;
        Params.repGrat         = 6;     % # of gratings stim repeats (total G's in V of NS stim protocol.
        
    case 'NS2'
        Params.durBlank        = 4;
        Params.durBlank_Grats  = 1;
        Params.durStim         = 13;
        Params.durGratings     = 1;
        Params.repGrat         = 6;     % # of gratings stim repeats (total G's in V of NS stim protocol.
        
end

%% Stim variables..................................................
TotalDurStim = data.stimDur;

NumSamplsNew = TotalDurStim.*data.newFrameRate;

TimeEachTrial   =  (4*17)+(8*12);
SamplsEachTrial = TimeEachTrial*data.newFrameRate;
NumTrials       = NumSamplsNew./SamplsEachTrial;

SamplsMovie     = 4*17*data.newFrameRate;
SamplsGratings  = 8*12*data.newFrameRate;

NumRepsMov  = 4;

%% Partition data..................................................
% Partitions into stimulus response frames
tic
ct = 0;
for i = 1:size(data.dFF,1)
    for j = 1:size(data.dFF,2)
        ct = ct+1;
        clc
        fprintf('Partitioning data...\n');
        fprintf('Completed:%0.3f%% (%0.2fs)\n',(ct/(size(data.dFF,2)^2))*100,toc);
        %         progressbar(ct/(size(data.dFF,2)^2));
        dFFx = squeeze(data.dFF(i,j,:) );
        Trials = reshape(dFFx, [SamplsEachTrial, NumTrials] );
        
        for tr = 1:NumTrials
            foo = squeeze( Trials(:,tr) );
            NatImages{tr}(i,j,:) =  foo(1:SamplsMovie);
            Gratings{tr}(i,j,:)  =  foo(SamplsMovie+1:SamplsMovie+SamplsGratings);
        end;
        
    end;
end;

% Save mat file containing partitioned dff
if saveFlag
    fprintf('Saving...NatImages...\n');
    fname_nat=[fname(1:end-8),'NatImages'];
    save(fname_nat,'NatImages');
    fprintf('Saving...NatImages...Done!\n');
    fprintf('Saving...Gratings...\n');
    fname_grat=[fname(1:end-8),'Gratings'];
    save(fname_grat,'Gratings', '-v7.3');
    fprintf('Saving...Gratings...Done!\n');
end

%% Display Movie (optional)........................................
if displayFlag
    for fr = 1:SamplsMovie
        figure(1); set(gcf,'color','w');
        subplot(3,2,1); imagesc( NatImages{1}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,2); imagesc( NatImages{2}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,3); imagesc( NatImages{3}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,4); imagesc( NatImages{4}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,5); imagesc( NatImages{5}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,6); imagesc( NatImages{6}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        drawnow;
    end;
    
    for fr = 1:SamplsGratings
        figure(2); set(gcf,'color','w');
        subplot(3,2,1); imagesc( Gratings{1}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,2); imagesc( Gratings{2}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,3); imagesc( Gratings{3}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,4); imagesc( Gratings{4}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,5); imagesc( Gratings{5}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        subplot(3,2,6); imagesc( Gratings{6}(:,:,fr) ); axis square; axis off; caxis([0,500]);
        drawnow;
    end;
end;

%% Compute NS Reliability map......................................
tic
fprintf('Computing Nat Reliability...\n');

DurMov = Params.durStim;
DurBl  = Params.durBlank;
P  = nchoosek( 1:24, 2);

ct = 0;
close all

for i = 1:size(data.dFF,1)
    for j = 1:size(data.dFF,2)
        ct = ct+1;
        clc
        fprintf('Computing Nat Reliability....\n');
        fprintf('Completed:%0.3f%% (%0.2f min)\n',(ct/(size(data.dFF,2)^2))*100,toc/60);
        S_trials = NaN((DurMov+DurBl)*data.newFrameRate,NumRepsMov);
        Sx = [];
        
        % make raster plot for pixel.......................................
        for tr = 1:NumTrials
            S            = squeeze(NatImages{tr}(i,j,:));
            S_trials     = (reshape( S', [(DurMov+DurBl)*data.newFrameRate,NumRepsMov] ))';
            Sx           = cat(1, Sx, S_trials);
        end;
        
        % evaluate Reliability for that pixel..............................
        for perm = 1:size(P,1);
            T1 = P(perm,1);
            T2 = P(perm,2);
            R(perm) =  distcorr( Sx(T1, DurBl*data.newFrameRate + 1: end)', Sx(T2,DurBl*data.newFrameRate + 1: end)' );
        end;
        RelNat(i,j) = nanmean(R);
        clear R Sx S S_trials;
    end;
end;
dur(1) = toc;
fprintf('...Done!\n Elapsed time: %0.2fmins\n',dur(1)/60);

fprintf('Saving...Nat Reliability...\n');
fname_rel=[fname(1:end-8),'RelNat'];
save(fname_rel,'RelNat');
fprintf('Saving...Nat Reliability...Done!\n');

%% Compute Grat Reliability map....................................
tic
fprintf('Computing Gratings Reliability...\n');

DurGrat = Params.durGratings;
DurBl  = Params.durBlank_Grats;
NumReps = Params.repGrat;
NumOris  = 12;
P       = nchoosek( 1:NumReps, 2);

ct = 0;
close all

for i = 1:size(data.dFF,1)
    for j = 1:size(data.dFF,2)
        ct = ct+1;
        clc
        fprintf('Computing Grat Reliability....\n');
        fprintf('Completed:%0.3f%% (%0.2f min)\n',(ct/(size(data.dFF,2)^2))*100,toc/60);
        S_trials = NaN((DurGrat+DurBl)*data.newFrameRate,NumOris);
        
        for ori = 1:12
            Sx{ori} = [];
        end;
        
        switch protocol
            case 'NS2'
            % make raster plot for pixel.......................................
            for tr = 1:Params.repGrat
                S            = squeeze(Gratings{tr}(i,j,:));
                S_trials     = (reshape( S', [(DurGrat+DurBl)*data.newFrameRate,NumOris*4] ))';
                for ori = 1:NumOris
                    idxOri = [ori:12:NumOris*4];
                    foo  =  squeeze(S_trials(idxOri,:));
                    Sx{ori} = cat(1, Sx{ori}, foo);
                    clear foo
                end;
            end;
            
            case 'NS1'
                % make raster plot for pixel.......................................
            for tr = 1:Params.repGrat
                S            = squeeze(Gratings{tr}(i,j,:));
                S_trials     = (reshape( S', [(DurGrat+DurBl)*data.newFrameRate,NumOris] ))';
                for ori = 1:NumOris
                    foo  =  squeeze(S_trials(ori,:));
                    Sx{ori} = cat(1, Sx{ori}, foo);
                    clear foo
                end;
            end;
                
            
        end; % end switch
        
%         Sxx = cell2mat(Sx');
%         
%         P = nchoosek(1:24,2);
%         for mc = 1:3
%             indx = randperm(72);
%             values = indx(1:24);
%             Sxx_this = Sxx(values,:);
%             for perm = 1:size(P,1);
%                  T1 = P(perm,1);
%                  T2 = P(perm,2);
%                  R(mc,perm) =  distcorr( Sxx_this(T1,DurBl*data.newFrameRate + 1: end)', Sxx_this(T2,DurBl*data.newFrameRate + 1: end)');
%             end;
%         end;
%         RelGrat(i,j) = mean(nanmean(R,2));
%         rawRelGrat = R;
%          clear R Sx S S_trials;   
        
        % evaluate Reliability for that pixel..............................
        for ori = 1:NumOris
            for perm = 1:size(P,1);
                T1 = P(perm,1);
                T2 = P(perm,2);
                R(ori,perm) =  distcorr( Sx{ori}(T1, DurBl*data.newFrameRate + 1: end)', Sx{ori}(T2,DurBl*data.newFrameRate + 1: end)' );
            end;
        end;
        RelGrat(i,j) = mean(nanmean(R,2)); %mean over perms, %max over oris.
        clear R Sx S S_trials;
    end;
end;

dur(2) = toc;
fprintf('...Done!\n Elapsed time: %0.2fmins\n',dur(2)/60);

fprintf('Saving...Grat Reliability...\n');
fname_rel=[fname(1:end-8),'RelGrat'];
save(fname_rel,'RelGrat');
% fname_rel=[fname(1:end-8),'rawRelGrat'];
% save(fname_rel,'rawRelGrat');
fprintf('Saving...Grat Reliability...Done!\n');
pause(1);

%% Plot Rel Maps...................................................

if plotFlag
    plotRelMaps(fname, RelNat, RelGrat);
end

%% Done!!!!!

fprintf('NS anlaysis duration: %0.2fmins\n',dur(1)/60);
fprintf('Grating analysis duration: %0.2fmins\n',dur(2)/60);
fprintf('Total elapsed time: %0.2fmins\n',sum(dur)/60);

