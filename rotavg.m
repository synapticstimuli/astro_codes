% rotavg.m - function to compute rotational average of (square) array
% 
% function f = rotavg(array)
%
% array can be of dimensions N x N x M, in which case f is of 
% dimension NxM.  N should be even.
%

function f = rotavg(array)

[N N M]=size(array);

[X Y]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);

[theta rho]=cart2pol(X,Y);

% fprintf('initializing indices\n');

rho=round(rho);
if mod(N,2)~= 0
    i=cell( (N+1)/2,1);
    f=zeros( (N+1)/2,M);
else
    i=cell(N/2+1,1);
    f=zeros(N/2+1,M);

end;
for r=0:N/2
  i{r+1}=find(rho==r);
end

% fprintf('doing rotational average\n');

for m=1:M

  a=array(:,:,m);
  for r=0:N/2
    f(r+1,m)=mean(a(i{r+1}));
  end
  
end