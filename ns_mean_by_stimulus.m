function dataOut = ns_mean_by_stimulus(NS, dataIn)

condition_num = ns_get(NS, 'condition_num');
conditions = unique(condition_num);
sz = size(dataIn); ndim = sum(sz>1); 
if ndim == 1, dataOut = NaN(length(conditions), 1);
else          dataOut = NaN(size(dataIn,1), length(conditions)); end


for ii = 1:length(conditions)
    if ndim == 1,
        these_data = dataIn(condition_num==conditions(ii),:);
        dataOut(ii) = mean(these_data(:));    
    else
        these_data = dataIn(:,condition_num==conditions(ii),:);
        these_data = reshape(these_data, sz(1), []);
        dataOut(:,ii) = mean(these_data,2);    
    end
    
end

return