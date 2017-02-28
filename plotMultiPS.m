function [handle] = plotMultiPS(maps, handle)
if nargin==0
    tempData = {natThreshImg;gratCCImg;sponTreshImg};
else
    tempData = maps;
end
%%
if iscell(tempData)
    numData = length(tempData);
else
    numData = 1;
    tempData = {maps};
end

if ~exist('handle','var')
    figure; hold on;
    handle = gca;
end
colors = [0 0 0; 1 0 0; 0 0 1];

for p = 1:numData
    currData = tempData{p};
    faxis = 0:size(currData,1)/2;
    
    AX = fftshift(fft2(currData));
    PS = abs(AX).^2;
    
    PsY = rotavg(PS);
    
    %         if exist('currHandle','var')
    %             figure(handle);
    %             figHandle = gcf;
    %         else
    %             figHandle = figure; set(gcf,'color','w');
    %         end
    loglog(faxis, PsY,'Color',colors(p,:),'linewidth',3); hold on;
    axis square; box off;
%     set(gca,'XLim', [2 11], 'YScale','log');
    set(gca,'XLim', [2 11], 'YScale','log');
    %         % The following will modify the x axis if provided a range during the
    %         % initial funciton call (eg, [2 11]). Otherwise, it will leave the XLim property as is
    %         % (full range based on data set).
    %         if axisFlag
    %             if isempty(axisMod)
    %                 axisMod = get(gca,'XLim');
    %                 set(gca,'XLim', axisMod);
    %             else
    %                 set(gca,'XLim', axisMod);
    %             end
    %         end
    
    % Change default YLims to be a minimu of 100 (10^2)
    ylims = get(gca,'YLim');
    if ylims(2) < 100
        set(gca,'YLim',[1 100]);
    else
        set(gca,'YLim',[1 100000]);
    end
    xlabel('Spatial Freq (cyc / pix)');
    ylabel('log_1_0 Power');
    title('Power Spectra')
    % Additional axes properties
    set(gca,'XMinorTick','off','YMinorTick','off','XTick',[1 11],...
        'FontSize',12,'XTickLabel',{'1';'10'});
end
end
