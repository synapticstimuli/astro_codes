function [logfreq, logPf, gradient] = plotPowerSpectrumRVR(x,figNum)
% Plot power spectrum of image/movies. Generates coefficient for curve
% (gradient) for each PS.

%% Input args
if nargin == 2
    figure(figNum); 
    hold on;
elseif nargin == 1
    figNum = figure;
    hold on;
end

% gray = [.7,.7,.7];
% cmap = b2r(0,1);
% linecolor=cmap(20,:);

linemap = lines(5);

if length(size(x))  == 3
    numFrames = size(x,3);
    numPx     = size(x,1);
else
    numFrames = 1;
    numPx     = size(x,1);
end;

% convFact =  0.72;
 convFact = 1;
f       = -numPx/2 : numPx/2-1;
freq    = [0 : numPx/2]'.*convFact;
logfreq = log10(freq(2:numPx/2));

for i = 1:numFrames
    im   =  squeeze( x(:,:,i) );
    imf  = fftshift( fft2(im, numPx, numPx) ) ;
    impf = abs(imf).^2;
    Pf  = rotavg(impf);
    logPf(:,i) = log10(Pf(2:numPx/2));
end;

for i = 1:size(x,3)
% perform best fit
[Coeff, Resnorm, residual] = lsqcurvefit( @objFun, [-1,0], logfreq(2:10), logPf(2:10,i) );
gradient(i) = Coeff(1);
end;
 
plot( logfreq(2:10), mean(logPf(2:10),2), 'linewidth',3); hold on;



