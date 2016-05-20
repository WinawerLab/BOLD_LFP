function [dataOut,dataStd] = ns_mean_by_stimulus(NS, dataIn,varargin)

% usage:
% [dataOut,dataStd] = ns_mean_by_stimulus(NS, ns_get(NS, 'alpha')),'all')

condition_num = ns_get(NS, 'condition_num');
conditions = unique(condition_num);
sz = size(dataIn); ndim = sum(sz>1); 
if ndim == 1, 
    dataOut = NaN(length(conditions), 1);
    dataStd = NaN(length(conditions), 1);
else
    dataOut = NaN(size(dataIn,1), length(conditions)); 
    dataStd = NaN(size(dataIn,1), length(conditions)); 
end

if ~isempty(varargin) % take all data, even or odd trials
    if isequal(varargin{1},'odd') % take odd trials
        for ii = 1:length(conditions)
            if ndim == 1,
                these_data = dataIn(condition_num==conditions(ii),:);
                dataOut(ii) = mean(these_data(1:2:end));
                dataStd(ii) = std(these_data(1:2:end));
            else
                these_data = dataIn(:,condition_num==conditions(ii),:);
                these_data = reshape(these_data, sz(1), []);
                dataOut(:,ii) = mean(these_data(1:2:end),2);
                dataStd(:,ii) = std(these_data(1:2:end),[],2);
            end
        end
    elseif isequal(varargin{1},'even') % take even trials
         for ii = 1:length(conditions)
            if ndim == 1,
                these_data = dataIn(condition_num==conditions(ii),:);
                dataOut(ii) = mean(these_data(2:2:end));
                dataStd(ii) = std(these_data(2:2:end));
            else
                these_data = dataIn(:,condition_num==conditions(ii),:);
                these_data = reshape(these_data, sz(1), []);
                dataOut(:,ii) = mean(these_data(2:2:end),2);
                dataStd(:,ii) = std(these_data(2:2:end),[],2);
            end
         end
    elseif isequal(varargin{1},'all') % take all trials
        for ii = 1:length(conditions)
            if ndim == 1,
                these_data = dataIn(condition_num==conditions(ii),:);
                dataOut(ii) = mean(these_data(:));
                dataStd(ii) = std(these_data(:));
            else
                these_data = dataIn(:,condition_num==conditions(ii),:);
                these_data = reshape(these_data, sz(1), []);
                dataOut(:,ii) = mean(these_data,2);
                dataStd(:,ii) = std(these_data,[],2);
            end
        end
    end
elseif isempty(varargin) 
    for ii = 1:length(conditions)
        if ndim == 1,
            these_data = dataIn(condition_num==conditions(ii),:);
            dataOut(ii) = mean(these_data(:));    
            dataStd(ii) = std(these_data(:));    
        else
            these_data = dataIn(:,condition_num==conditions(ii),:);
            these_data = reshape(these_data, sz(1), []);
            dataOut(:,ii) = mean(these_data,2);    
            dataStd(:,ii) = std(these_data,[],2);    
        end

    end
end
return