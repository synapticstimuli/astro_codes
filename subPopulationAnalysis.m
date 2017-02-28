function figHandles = subPopulationAnalysis(gradientNat, gradientGrat,natRelData,gratRelData)
if nargin==0
    load ('PSdata.mat');
    load('NatRelData.mat');
    load('GratRelData.mat');
end
    
%% Parse data

supMGratIdx = gradientGrat > mean(gradientGrat);
supMNatIdx = gradientNat > mean(gradientNat);

subMGratIdx = gradientGrat < mean(gradientGrat);
subMNatIdx = gradientNat < mean(gradientNat);

%% Sup mean data
supFiles = {};
idx = find(supMNatIdx);
for i=1:length(idx)
    cur = idx(i);
    supFiles{i,1} = natRelData{cur}.fname;
    supdata(i,1) = abs(gradientNat(cur));
end

idx = find(supMGratIdx);
for i=1:length(idx)
    cur = idx(i);
    supFiles{i,2} = gratRelData{cur}.fname;
    supdata(i,2) = abs(gradientGrat(cur));
end

%% Sub mean data
subFiles = {};
idx = find(subMNatIdx);
for i=1:length(idx)
    cur = idx(i);
    subFiles{i,1} = natRelData{cur}.fname;
    subdata(i,1) = abs(gradientNat(cur));
end

idx = find(subMNatIdx);
for i=1:length(idx)
    cur = idx(i);
    subFiles{i,2} = gratRelData{cur}.fname;
    subdata(i,2) = abs(gradientGrat(cur));
end

%% Check files are the same for grat and nats

for i = 1:length(subFiles)
    gfile = subFiles{i,1};
    nfile = subFiles{i,2};
    fileChk(i) = strcmp(gfile,nfile);
end

if ~sum(fileChk)
    fprintf('All subM files the same.\n');
end

for i = 1:length(supFiles)
    gfile = supFiles{i,1};
    nfile = supFiles{i,2};
    fileChk(i) = strcmp(gfile,nfile);
end
if ~sum(fileChk)
    fprintf('All supM files the same.\n');
end

%% Plot sub data
cmap = [.3 .3 .3; 1 .3 .3; 0 0 0; 1 0 0];

f(1) = figure; hold on;
subplot(2,2,[1 2]); hold on
h(1) = histogram(subdata(:,1),10,...
    'FaceColor', 'k');        %Sub m Nat values
h(2) = histogram(subdata(:,2),10,...
    'FaceColor', 'r');        %Sub m Grat values
h(5) = histogram(supdata(:,1),10,...
    'FaceColor', 'k');        %Supra m Nat values
h(6) = histogram(supdata(:,2),10,...
    'FaceColor', 'r');        %Supra m Nat values

yline = get(gca,'ylim');
natXline = repmat(mean(abs(gradientNat)),1,2);
h(3) = line(natXline,yline,'LineWidth',6,'LineStyle',':','Color',cmap(1,:));

gratXline = repmat(mean(abs(gradientGrat)),1,2);
h(4) = line(gratXline,yline,'LineWidth',6,'LineStyle',':','Color',cmap(2,:));

set(h(1:2), 'FaceAlpha', 0.3)
set(h(5:6), 'FaceAlpha', 0.7)
dataNames = {'Nat', 'Grat','mean Nat','mean Grat'};
legend(dataNames,'Location','northwest');

% Plot boxScatter
subplot(2,2,3); hold on;
boxScatter(subdata,f(1));
title('Supra-mean population');
subplot(2,2,4)
boxScatter(supdata,f(1));
title('Sub-mean population');

set(f, 'Position', [200 50 940 840]);
%%
figHandles = f;
        