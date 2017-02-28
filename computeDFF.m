function final_dFF = computeDFF(F, sampling_rate)

type = 'Method1';

switch type
    case 'Method1'
        if nargin<2; sampling_rate = 6; end;    
        
        % convert f to dff
        Fb = prctile(F(:),25);
        if Fb
            dFF = (F-Fb)./Fb;
        else
            dFF = F * 0;
        end
        
        detrend_dFF = (dFF);
        
        % filter
        if sampling_rate > 5
            span = ceil(sampling_rate / 2);   % Size of the averaging window
            window = ones(span,1)/span;
            smoothed_dFF = filter(window,1,detrend_dFF);
        else
            smoothed_dFF = detrend_dFF;
        end
        
        % bring back to 0
        min_F = prctile(smoothed_dFF(:),10);
        final_dFF = smoothed_dFF - min_F;
        
        % scale to percentage
        final_dFF = final_dFF * 100;
        
    case 'Method2'
        [KSD,Xi] = ksdensity(F);
        [~,maxIdx]= max(KSD);
        F0 = Xi(maxIdx);
        final_dFF =  ((F-F0)/F0);
        
         detrend_dFF = (final_dFF);
         % filter
        if sampling_rate > 5
            span = ceil(sampling_rate / 2);   % Size of the averaging window
            window = ones(span,1)/span;
            smoothed_dFF = filter(window,1,detrend_dFF);
        else
            smoothed_dFF = detrend_dFF;
        end
            % bring back to 0
        min_F = prctile(smoothed_dFF(:),10);
        final_dFF = smoothed_dFF - min_F;
        
        % scale to percentage
        final_dFF = final_dFF * 100;
end