function rm = RiskMetrics(data,lamda,method);
%---------------------------------------------------
% PURPOSE: 
% Estimates the univariate and multivariate case of RiskMetrics
%---------------------------------------------------
% USAGE: 
% rm = RiskMetrics(data,alpha)
%
% INPUTS: 
% data = ( m x n ) vector
% lamba = the scale parameter
% method = Univariate or Multivariate
% 
% OUTPUTS:
% rm = ( m x n ) volatility vector for the univariate case or 
% an [( n x n )x m] covariance vector for the multivariate case
% 
% NOTES:
% In case no lamda parameter is defined, the 0.95 is selected
%---------------------------------------------------
% Author:
% Alexandros Gabrielsen, a.gabrielsen@city.ac.uk
% Date: 09/2010
%---------------------------------------------------
if nargin == 0 
    error('Data, Lamda, Method') 
end

if isempty(lamda) 
    lamda = 0.95;
end

[a b] = size(data);
rm = [];

if isequal(method, 'Univariate')
   rm(1:b) = std(data); % starting point
for i = 1 : b 
    for j = 2 : a 
        rm(j,i) = sqrt(lamda*rm(j-1,i)^2 + (1-lamda)*data(j-1,i)^2);      
    end  
end
elseif isequal(method, 'Multivariate')
    % Engle and Shepphard (2001) M-RiskMetrics
    rm(:,:,1) = cov(data); % starting point
    for ii = 2 : a
    rm(:,:,ii) = (1-lamda)*data(ii-1,:)'*data(ii-1,:) + (lamda)*rm(:,:,ii-1);
    end
else disp('Select Method: Univariate or Multivariate');
end
end