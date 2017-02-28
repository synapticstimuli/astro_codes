function [rois,NatImages,RelNat,Gratings]=vizReliabilityMap_RG2(fname)

% Processes the data file and calculates the reliability of the responses
% to both ns and gratings. Distuinguishes between different stim protocols 
% (NS1 or NS2). Runs about 1.5 hrs per file. From RVR.

%% Load data file
clc
if nargin==0
    fprintf('Select data.mat file...\n')
    [fname, pathname] = uigetfile('*data.mat','Load .mat file containing data structure...');
    cd(pathname);
    fprintf('Loading file...\n');
    load (fname, 'data')
%     displayData(data);
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

%% Read protocol from file name

% User variables
displayFlag     =   0;
saveFlag        =   1;
plotFlag        =   0;
redoAnalysis    =   0;
nsProtocol      =   'NS1';       %Select between original NS or NS2.
stimProtocol    =   'Astro';      % Select between Astro or Soma
skipNM          =   0;
skipGrat        =   1;

%% Stimulus parameters
Params.repTrials = 6;                   % Num of trials (movies + gratings)
Params.repMovies = 4;                   % Num of movie reps per trial
Params.numOris = 12;                    % Num og orientations per trial

switch nsProtocol
    case 'NS1'
        switch stimProtocol
            case 'Astro'
                data.stimDur           = 984;
                Params.durBlank        = 4;
                Params.durBlank_Grats  = 6;
                Params.durStim         = 13;
                Params.durGratings     = 2;
            case 'Soma'
                data.stimDur           = 792;
                Params.durBlank        = 2;
                Params.durBlank_Grats  = 4;
                Params.durStim         = 13;
                Params.durGratings     = 2;
        end
    case 'NS2'
        data.stimDur           = 984;
        Params.durBlank        = 4;
        Params.durBlank_Grats  = 1;
        Params.durStim         = 13;
        Params.durGratings     = 1;
end
fprintf('Protocol: %s\n', nsProtocol);
durTotal = Params.repTrials*(((Params.durBlank+Params.durStim )*Params.repMovies)...
           +((Params.durBlank_Grats+Params.durGratings)*Params.numOris));
fprintf('Stim = %s\nTotal Time = %d\n', stimProtocol,durTotal);
pause(3);

%% HouseKeeping
TotalDurStim = data.stimDur;

NumSamplsNew = TotalDurStim.*data.newFrameRate;

TimeEachTrial   = (Params.repMovies*(Params.durStim+Params.durBlank)) + (Params.numOris*(Params.durGratings+Params.durBlank_Grats));
SamplsEachTrial = TimeEachTrial*data.newFrameRate;
NumTrials       = NumSamplsNew./SamplsEachTrial;

SamplsMovie     = (Params.repMovies*(Params.durStim+Params.durBlank))*data.newFrameRate;
SamplsGratings  = (Params.numOris*(Params.durGratings+Params.durBlank_Grats))*data.newFrameRate;

%% Partition data
tic
ct = 0;
for i = 1:size(data.dFF,1)
    for j = 1:size(data.dFF,2)
        ct = ct+1;
        clc
        fprintf('Partitioning data...\n');
        fprintf('Completed:%0.3f%% (%0.2fs)\n',(ct/(size(data.dFF,2)^2))*100,toc);
        dFFx = squeeze(data.dFF(i,j,:) );
        Trials = reshape(dFFx, [SamplsEachTrial, NumTrials] );
        
        for tr = 1:NumTrials
            foo = squeeze( Trials(:,tr) );
            NatImages{tr}(i,j,:) =  foo(1:SamplsMovie);
            Gratings{tr}(i,j,:)  =  foo(SamplsMovie+1:SamplsMovie+SamplsGratings);
        end;
        
    end;
end;

pause(2)

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

%% Compute NS Reliability map
if ~skipNM
    tic
    fprintf('Computing Nat Reliability...\n');
    
    DurMov = Params.durStim;
    DurBl  = Params.durBlank;
    P  = nchoosek( 1:(Params.repTrials*Params.repMovies), 2);
    
    ct = 0;
    close all
    
    for i = 1:size(data.dFF,1)
        for j = 1:size(data.dFF,2)
            ct = ct+1;
            clc
            fprintf('Computing Nat Reliability....\n');
            fprintf('Completed:%0.3f%% (%0.2f min)\n',(ct/(size(data.dFF,2)^2))*100,toc/60);
            S_trials = NaN((DurMov+DurBl)*data.newFrameRate,Params.repMovies);
            Sx = [];
            
            % make raster plot for pixel.......................................
            for tr = 1:NumTrials
                S            = squeeze(NatImages{tr}(i,j,:));
                S_trials     = (reshape( S', [(DurMov+DurBl)*data.newFrameRate,Params.repMovies] ))';
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
else
    fprintf('NM rel analysis skipped...\n');
end

%% Compute Grat Reliability map
if ~skipGrat
    tic
    fprintf('Computing Gratings Reliability...\n');
    
    DurGrat = Params.durGratings;
    DurBl  = Params.durBlank_Grats;
    P       = nchoosek( 1:Params.repTrials, 2);
    
    ct = 0;
    close all
    
    for i = 1:size(data.dFF,1)
        for j = 1:size(data.dFF,2)
            ct = ct+1;
            clc
            fprintf('Computing Grat Reliability....\n');
            fprintf('Completed:%0.3f%% (%0.2f min)\n',(ct/(size(data.dFF,2)^2))*100,toc/60);
            S_trials = NaN((DurGrat+DurBl)*data.newFrameRate,Params.numOris);
            
            for ori = 1:12
                Sx{ori} = [];
            end;
            
            switch nsProtocol
                case 'NS2'
                    % make raster plot for pixel.......................................
                    for tr = 1:Params.repTrials
                        S            = squeeze(Gratings{tr}(i,j,:));
                        S_trials     = (reshape( S', [(DurGrat+DurBl)*data.newFrameRate,Params.numOris*4] ))';
                        for ori = 1:Params.numOris
                            idxOri = [ori:12:Params.numOris*4];
                            foo  =  squeeze(S_trials(idxOri,:));
                            Sx{ori} = cat(1, Sx{ori}, foo);
                            clear foo
                        end;
                    end;
                    
                case 'NS1'
                    % make raster plot for pixel.......................................
                    for tr = 1:Params.repTrials
                        S            = squeeze(Gratings{tr}(i,j,:));
                        S_trials     = (reshape( S', [(DurGrat+DurBl)*data.newFrameRate,Params.numOris] ))';
                        for ori = 1:Params.numOris
                            foo  =  squeeze(S_trials(ori,:));
                            Sx{ori} = cat(1, Sx{ori}, foo);
                            clear foo
                        end;
                    end;
                    
            end; % end switch
            
            % evaluate Reliability for that pixel..............................
            for ori = 1:Params.numOris
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
    fprintf('Saving...Grat Reliability...Done!\n');
    pause(1);
else
    fprintf('Grat rel  skipped...\n');
end

%% Plot Rel Map

if plotFlag
    plotVizRelMapPS_RG(fname, RelNat, RelGrat);
end
% %% Finish
% 
% fprintf('NS anlaysis duration: %0.2fmins\n',dur(1)/60);
% fprintf('Grating analysis duration: %0.2fmins\n',dur(2)/60);
% fprintf('Total elapsed time: %0.2fmins\n',sum(dur)/60);

