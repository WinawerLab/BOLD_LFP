function r2 = ns_cod(pred, data, varargin)

% pred = pred(:);
% data = data(:);

if ~isempty(varargin) 
    normalize = varargin{1}; 
else
    normalize = false;
end

if normalize
    pred = pred - mean(pred);
    data = data - mean(data);
    pred = pred / sum(pred.^2);
    data = data / sum(data.^2);
end
    
SSres = sum((pred-data).^2);
SSdat = sum((data-mean(data)).^2);

r2 = 1 - SSres/SSdat;
