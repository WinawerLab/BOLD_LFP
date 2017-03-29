function data_epoch = ecog_regresserp(data_epoch,t_base,stims)
% function to regress the ERP from some data
% USAGE:
% data_epoch = ecog_regresserp(data_epoch,t_base,stims)
% 
% data_epoch % epochs X time
% t_base % indices of the baseline period
% stims % different code for different conditions

% baseline correct
for m=1:size(data_epoch,1)%epochs
    x=squeeze(data_epoch(m,:));
    x=x-mean(x(t_base));
    data_epoch(m,:)=x;
end


% regress erp out per condition
for m=1:size(data_epoch,1)%epochs
    x=data_epoch(m,:);
    % regress ERP out
    s=stims(m);
    av_erp=squeeze(mean(data_epoch(stims==s,:),1));
    [~,~,reg_R] = regress(x',av_erp');
    data_epoch(m,:)=reg_R;
end
