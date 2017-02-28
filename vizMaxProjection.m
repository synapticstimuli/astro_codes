function [MaxProj] = vizMaxProjection(dFF)

% fprintf('Computing the max projection...');
% for i = 1:size(dFF,1)
%     for j = 1:size(dFF,1)
%         MaxProj = max( squeeze( dFF(i,j,:)));
%     end;
% end;
MaxProj = max(dFF, [], 3);
fprintf('...Done\n');