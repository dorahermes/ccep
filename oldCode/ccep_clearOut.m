function data = ccep_clearOut(data,outlier_matrix)
% Sets outliers in the data to NaN
% clearOut(data,outlier_matrix);
%
% Inputs:
%   - data is epoched ccep data (elXtimeXepoch)
%   - outlier matrix has a 1 for each outlier, 0 for none (elXepoch)
%
%
% Hermes and Vassileva, 2016, UMC Utrecht

for el=1:size(data,1)
    for l=1:size(data,3)
        if outlier_matrix(el,l)==0
            continue
        else
            data(el,:,l) = NaN; 
        end
    end
    disp (['el ' int2str(el) ' - epochs with outliers cleared.'])
end