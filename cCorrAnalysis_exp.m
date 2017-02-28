function cCorrAnalysis_exp(stimType)
progressbar;
%% Input params
rVal = 0.5;
ccFolder = 'D:\Users\Rodrigo\TempData\aTetO\ccMap_Data';
plotPS = 0;
plotHist = 1;
saveFlag = 1;
%% Load files
% Load vis response data
switch stimType
    case 'NS-'
        try
%             pathname = [pwd filesep];
            [foo,curDir,foo] = fileparts(pwd);
            fname = dir('*NatImages.mat');
            if isempty(fname)
                fprintf('NatImages.mat file not found in %s...\n',curDir);
                pause(1);
                return
            end
            natName =fname.name;
            gratName = strrep(natName,'NatImages','Gratings');
            fprintf('Loading %s:...\n', natName);
            load(natName);
            fprintf('Loading %s:...\n', gratName);
            load(gratName);
            expName = curDir;
            
        catch
            fprintf('Select NatImages data.mat file...\n\n')
            [fname, pathname] = uigetfile('*NatImages.mat','Load NatImages.ma file...');
            cd(pathname);
            fprintf('Loading %s:...\n', fname(1:end-23));
            fprintf('NS data...\n');
            load([pathname fname]);
            fprintf('Gratings data...\n');
            load([pathname fname(1:end-13), 'Gratings']);
        end
    case 'Spo-'
        try
            [foo,curDir,foo] = fileparts(pwd);
            fname = dir('*data.mat');
            fnameSpon = fname.name;
            fprintf('Loading %s:...\n', fnameSpon);
            load(fnameSpon);
            
            expName = [curDir '-Spo'];
            sponDff = data.dFF;
        catch
            % Load spontaneous data
            fprintf('Select spontaneous data.mat file...\n\n')
            [fnameSpon, pathnameSpon] = uigetfile('*data.mat','Load data.mat file containing spontaneous data...');
            cd(pathnameSpon);
            fprintf('Spontaneous data...\n');
            load (fnameSpon, 'data')
            expName = fnameSpon(1:end-18);
            sponDff = data.dFF;
        end
end

%% Extract dff data    
switch stimType
    case 'NS-'
        % NS data
        natdff = [];
        for i =1:6
            natdff = cat(3,natdff,NatImages{i});
        end
        natdff = cropMat(natdff,2);
        
        % Grat data
        gratdff = [];
        for i =1:6
            gratdff = cat(3,gratdff,Gratings{i});
        end
        gratdff = cropMat(gratdff,2);
        
    case 'Spo-'
        % Spontaneous data
        sponDff = cropMat(sponDff,2);
end

%% Run correlation and thresholding
switch stimType
    case 'NS-'
        [natCCImg,natThreshImg] = CrossCorrImage(rVal, natdff);
        set(gcf, 'Name',[expName '-NatCC']);
        [gratCCImg,gratThreshImg] = CrossCorrImage(rVal, gratdff);
        set(gcf, 'Name',[expName '-GratCC']);
    case 'Spo-'
        [sponCCimg, sponTreshImg] = CrossCorrImage(rVal, sponDff);
        set(gcf, 'Name',[expName 'CC']);
end
%% Plot PS of cc maps 
if plotPS
    % h1 = figure; hold on
    % plotMultiPS({gratThreshImg,natThreshImg,sponTreshImg});
    dataSets.mats = {gratThreshImg,natThreshImg,sponTreshImg};
    dataSets.titles = {'Grat','NM','Spon'};
    for i = 1:length(dataSets.mats)
        plotMultiPS(dataSets.mats{i})
        title(['PS - ' dataSets.titles{i}]);
        legend('off');
    end
end

%% Histogram of pixels in cc maps
if plotHist
    switch stimType
        case 'NS-'
            f(1)    = figure('visible','off','name',expName);
            hst(1)  = subplot(1,2,1);
            h(1)    = histogram(natCCImg(natCCImg>=0.02));
            title('Nat Movies');
            hst(2)  = subplot(1,2,2);
            h(2)    = histogram(natCCImg(natCCImg>=rVal));
            title(['Nat Movies, R>(', num2str(rVal),')']);
            suptitle(expName);
            
            % Set figure properties
            ylims = max(cell2mat(get(hst,'YLim')));
            xlims = max(cell2mat(get(hst,'XLim')));
            for i=1:length(hst);
                set(hst(i),'YLim',ylims,'XLim',xlims);
            end
            set(f,...
                'color'         ,'w'                ,...
                'Position'      , [101 544 800 270]) ;
            
            f(2)    = figure('visible','off','name',expName);
            hst(1)  = subplot(1,2,1);
            h(1)    = histogram(gratCCImg(gratCCImg>=0.02));
            title('Gratings');
            hst(2)  = subplot(1,2,2);
            h(2)    = histogram(gratCCImg(gratCCImg>=rVal));
            title(['Gratings, R>(', num2str(rVal),')']);
            suptitle(expName);
            
            % Set figure properties
            ylims = max(cell2mat(get(hst,'YLim')));
            xlims = max(cell2mat(get(hst,'XLim')));
            for i=1:length(hst);
                set(hst(i),'YLim',ylims,'XLim',xlims);
            end
            set(f,...
                'color'         ,'w'                ,...
                'Position'      , [101 544 800 270]) ;
            
        case 'Spo-'
            f(1)    = figure('visible','off','name',expName);
            hst(1)  = subplot(1,2,1);
            h(1)    = histogram(sponCCimg(sponCCimg>=0.02));
            title('Spontaneous');
            hst(2)  = subplot(1,2,2);
            h(2)    = histogram(sponCCimg(sponCCimg>=rVal));
            title(['Spontaneous, R>(', num2str(rVal),')']);
            suptitle(expName);
            
            % Set figure properties
            ylims = max(cell2mat(get(hst,'YLim')));
            xlims = max(cell2mat(get(hst,'XLim')));
            for i=1:length(hst);
                set(hst(i),'YLim',ylims,'XLim',xlims);
            end
            set(f,...
                'color'         ,'w'                ,...
                'Position'      , [101 544 800 270]) ;
    end
    
    
   
    
end


%% Save
ccMaps.name = expName;
switch stimType
    case 'NS-'
        ccMaps.natCC = natCCImg;
        ccMaps.gratCC = gratCCImg;
        ccMaps.natThr = natThreshImg;
        ccMaps.gratThr = gratThreshImg;
    case 'Spo-'
        ccMaps.sponCC = sponCCimg;
        ccMaps.sponThr = sponTreshImg;
end
if saveFlag
    fprintf('Saving data...\n');
    cd(ccFolder);
    save([expName '.mat'],'ccMaps');
    fprintf('Saving figures...\n');
%     t = get(0,'children');
    for i=1:length(f)
        saveas(f(i),[f(i).Name '.jpg'],'jpeg')
        saveas(f(i),[f(i).Name '.eps'],'eps2c')
    end
end
% pause(3)
% close all
% 
% cd(pathname);
% cd(['..' filesep])
% fprintf('Finished with exp: %s:\n', fname(1:end-23));
end