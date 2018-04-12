% myregress.m
%
% purpose: regression statistics
%   usage: [d] = myregress(x, y, m1,dispfit)
%          x = column vector of independent variable
%          y = column vector of dependent variable
%          m1 = slope to check difference of
%          m = slope of regression line
%          b = y-intercept
%          r2 = coefficient of determination
%          pm = p-value for slope
%          pb = p-value for y-intercept
%      by:  justin gardner
%    date:  9/3/98
function [d b r2 pm pb] = myregress(x, y,m1,dispfit)

if ((nargin <2) || (nargin > 4))
  help myregress, return,
end
if (nargin <3)
  m1 = 0;
end
if (nargin <4)
  dispfit = 0;
end
if (dispfit == 1)
  dispfit = 'k-';
end

% make into column vectors
if (size(x,1)==1),x=x';,end
if (size(y,1)==1),y=y';,end
% denan
goodvals = find(~isnan(x) & ~isnan(y));
x = x(goodvals);
y = y(goodvals);

n = length(x);

if (n <= 2) 
  d.b = nan;d.m = nan;d.r2 = nan;d.pm = nan;d.pb = nan;d.r = nan;
  return
end
% precalculate sum terms
sumx = sum(x);
sumy = sum(y);
sumx2 = sum(x.^2);
sumy2 = sum(y.^2);
sumxy = sum(x.*y);

% use formulas to calculate slope and y-intercept
m = (n*sumxy-sumx*sumy)/(n*sumx2-sumx^2);
b = (sumy*sumx2 - sumx*sumxy)/(n*sumx2 - sumx^2);

% find least squares fit for a line of the form y = m*x + b
%A = x;
%A(:,2) = 1;
%coef = ((A'*A)^-1)*A'*y;
%m = coef(1);
%b = coef(2);

% calculate r^2
num = ((m*x + b) - y).^2;
denom = (y-mean(y)).^2;
r2 = 1 - (sum(num)/sum(denom));

% calculate standard error of the estimate
Sx = std(x);
Sy = std(y);
Syx = sqrt(((n-1)/(n-2))*(Sy^2-(m^2)*Sx^2));

% calculate standard error of the slope
Sm = (1/sqrt(n-1))*(Syx/Sx);

% calculate standard error of the y-intercept
Sb = Syx * sqrt((1/n) + (mean(x)^2)/((n-1)*Sx^2));

if (n <= 2), pm = -1;, pb = -1;, return;, end;

% calculate t-statistic and p-value for slope
% note, use incomplete beta function because
% matlab's tpdf function gives strange results.
Tm = (m-m1)/Sm;
if isnan(Tm)
  pm = nan;
else
  pm = betainc((n-2)/((n-2)+Tm^2),(n-2)/2,.5);
end

% calculate t-statistic and p-value for y-intercept
Tb = b/Sb;
pb = betainc((n-2)/((n-2)+Tb^2),(n-2)/2,.5);

% get all stats into structure
d.m = m;d.b = b;d.r2 = r2;d.pm = pm;d.pb = pb;

% Straight-forward of computing regression coefficient
%d.r = 0;
%for i = 1:n
%  d.r = d.r+(x(i)-mean(x))*(y(i)-mean(y))
%end
%d.r = d.r/((n-1)*Sx*Sy);
% r is related to slope by:
d.r = d.m*Sx/Sy;

% plot the results on current axis
if (dispfit)
  a = axis;
  hold on
  plot([a(1) a(2)],b+m*[a(1) a(2)],dispfit);
end

if nargout > 1
  d = d.m;
end
