% function filterROIs(~)

[settings, button] = settingsdlg(...
    'Description'                                       ,'Set parameters:'                      ,...
    {'Response p value threshold','p'}                  ,'0.05'                                 ,...
    {'Threshold NM responses?','filtNMflag'}            ,[false]                                ,...
    {'Threshold Gratings responses?','filtGratflag'}    ,[false]                                ,...
    {'Threshold OSI responses?','filtOSIflag'}          ,[false]                                ,...
     {'OSI threshold','osiThresh'}                      ,'0.1'                                  ,...
    {'Threshold GOF?','filtGOFflag'}                    ,[false]                                ,...
    {'GOF threshold','gofThresh'}                       ,'0.5'                                  ...
    );

if strcmp(button,'cancel') || isempty(button)
    fprintf('\n*****************  No settings were entered. Script stopped.  *********************\n\n\n');
    return
end

%%

 NumROIs = numel(rois.Label);
 
 % filter rois based on the individual threshholds...
 sigROIs = [];
 sigROIs.NmResp = find(rois.sigNMResp < settings.p);
 sigROIs.GratResp = find(rois.sigOriResp < settings.p);
 sigROIs.OSI = find(rois.OSI > settings.osiThresh);
 sigROIs.GOF =find(rois.GOF > settings.gofThresh);
 % Extract the filter names
 sigNames = fieldnames(sigROIs);                                
 % Identify the filters to use
 filtIdx = [settings.filtNMflag, settings.filtGratflag,...      
     settings.filtOSIflag, settings.filtGOFflag];
 filtFlags = sigNames(find(filtIdx));
 
 filtROI = rois.Label;
 filtStr = [];
 for i = 1:length(filtFlags)
     filtROI = intersect(filtROI,eval(['sigROIs.' filtFlags{i}]));
     filtStr = [filtStr, sprintf('%s ', filtFlags{i})];
 end
 filtStr = strrep(filtStr(1:end-1),' ', ', ');
 numFilt = length(filtROI);
 
 %% Create table............................................................
numRows = length(filtROI);
colNames = {'Rois','Reliability','sigNMResp','OSI','PO','DSI','sigOriResp','GOF',};
T = [];
T = table(NaN);
T = repmat(T, numRows, 8);
T.Properties.VariableNames = colNames;

for k = 1:numRows
    T.Rois(k,1)          = rois.Label(filtROI(k));
    T.Reliability(k,1)   = rois.Reliability(filtROI(k));
    T.sigNMResp(k,1)     = rois.sigNMResp(filtROI(k));
    T.OSI(k,1)           = rois.OSI(filtROI(k));
    T.PO(k,1)            = rois.PO(filtROI(k));
    T.DSI(k,1)           = rois.DSI(filtROI(k));
    T.sigOriResp(k,1)    = rois.sigOriResp(filtROI(k));
    T.GOF(k,1)           = rois.GOF(filtROI(k));  
end

T.Properties.UserData = ['Filters: ', filtStr];
clc;
disp(T.Properties.UserData);
disp(T);
tableName = [rois.fname(1:end-8), '-filtTable.txt'];
% writetable(T, tableName);
disp(sigROIs);

%% Append filtered rois...
answer = input('Append results (y or n)?','s');
if strcmp(answer,'y')
    rois.filtROI = filtROI;
    fprintf('Saving results...\n');
    save(rois.fname, 'rois', '-append');
end

% end
