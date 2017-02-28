function [figHandle] = plotPowerSpectra(RelGrat, RelNat, axisMod, handle,fname)
% From Rajeev, modified to included option to modify axis of final plot.
% Third input argument should be XLim range (eg [2 11])

saveFlag = 0;

% Load files 
if nargin == 0
    filelist = dir;
    for i = 3: length(filelist)
        if ~isempty (strfind(filelist(i).name,'RelNat.mat'))
            load(filelist(i).name);
        end
        if ~isempty (strfind(filelist(i).name,'RelGrat.mat'))
            load(filelist(i).name);
        end
    end
    axisFlag = 1;
elseif nargin == 1
    filelist = dir;
    axisFlag = 1;
    axisMod = RelGrat;
    saveFlag = 1;
    for i = 3: length(filelist)
        if ~isempty (strfind(filelist(i).name,'RelNat.mat'))
            load(filelist(i).name);
            saveName = strrep(filelist(i).name,'RelNat.mat','PS');
        end
        if ~isempty (strfind(filelist(i).name,'RelGrat.mat'))
            load(filelist(i).name);
        end
    end
elseif nargin == 2
    axisFlag = 0;
elseif nargin == 3
    axisFlag = 1;
elseif nargin ==5;
    axisFlag = 1;
    currHandle = 1;
%     saveFlag = 1;
     saveName = strrep(fname,'data.mat','PS');
else
    axisFlag = 0;
end

%% Plot PS

faxis = 0:size(RelGrat,1)/2;

AX = fftshift(fft2(RelGrat));
PSGrat = abs(AX).^2;

AX = fftshift(fft2(RelNat));
PSNat = abs(AX).^2;

gratY = rotavg(PSGrat);
natY = rotavg(PSNat);

if exist('currHandle','var')
    figure(handle);
    figHandle = gcf;
else
    figHandle = figure; set(gcf,'color','w');
end
loglog(faxis, gratY,'r','linewidth',3); hold on;
loglog(faxis, natY,'k','linewidth',3); hold on;
axis square; box off;

% The following will modify the x axis if provided a range during the
% initial funciton call (eg, [2 11]). Otherwise, it will leave the XLim property as is
% (full range based on data set).
if axisFlag
    if isempty(axisMod)
        axisMod = get(gca,'XLim');
        set(gca,'XLim', axisMod);
    else
        set(gca,'XLim', axisMod);
    end
end

% Change default YLims to be a minimu of 100 (10^2)
ylims = get(gca,'YLim');
if ylims(2) < 100
    set(gca,'YLim',[0 100]);
end

% Additional axes properties
set(gca,'XMinorTick','off','YMinorTick','off','XTick',[1 11],...
    'FontSize',12,'XTickLabel',{'1';'10'});
legend('NS','Grat');%,'Location','southoutside','Orientation','horizontal');
legend('boxoff');

%% Save figure
if saveFlag
    title(num2str(saveName(1:end-12)),'FontSize',10,'FontWeight','normal');
    fprintf('Saving figure (eps and jpg)...\n')
    saveas(figHandle,[saveName,'.eps'], 'eps2c');
    saveas(figHandle,[saveName,'.jpg']);
end
pause(4);
% close 
