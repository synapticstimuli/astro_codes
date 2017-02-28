xdff = computeDFF(x, 7);
calTrace = x;
lentrace = length(calTrace);
slWindow = 100;
newDFF = zeros(1,lentrace);
ct = 0;
for i = 1:slWindow:lentrace-slWindow
%     temp_median = median(calTrace(i:i+slWindow));
%     newDFF(i:i+slWindow) = (calTrace(i:i+slWindow) - temp_median) ./ temp_median;
    [KSD,Xi] = ksdensity(calTrace(i:i+slWindow));
    [~,maxIdx]= max(KSD);
    F0 = Xi(maxIdx);
    newDFF(i:i+slWindow) =  ((calTrace(i:i+slWindow)-F0)/F0)*100;
end

newDFF(1,end-slWindow) = (calTrace(end-slWindow) - temp_median) ./ temp_median;
% newDFF = newDFF*100;

figure('Position', [3 433 1876 545]);
subplot(311)
plot(x, 'k');
subplot(312)
plot(xdff,'b');
subplot(313)
plot(newDFF, 'r');