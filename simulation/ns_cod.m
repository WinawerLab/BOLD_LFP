function cod = ns_cod(x,y,varargin)
% calculates coefficient of determination (COD)
% COD: 1 - SSres/SSdata

% DH 2017

if ~isempty(varargin)
    rescale_opt = varargin{1};
else
    rescale_opt = 0;
end

% option to subtract mean and vector length normalize data, do this only if
% the data and estimate are not in the same units
if rescale_opt    
    % subtract mean
    x = x-mean(x);
    y = y-mean(y);

    % vector length normalize:
    x = x/sqrt(sum(x.^2));
    y = y/sqrt(sum(y.^2));
end

% here calculate 

SSres = sum((x - y).^2);

SSdata = sum((y - mean(y)).^2);

cod = 1 - (SSres/SSdata);