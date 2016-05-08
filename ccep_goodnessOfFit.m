%% Check if the model fits the data well.
% Checks relative root mean square error, Spearman's rho and p-value for
% Spearman's.
% Requires fitting parameters - ccep_setFitParms(t) and anatomical
% locations - ccep_anat_locs
% 
% Outputs:
% goodness_mat - el by elm by 3 matrix
% el - stimulated electrode
% elm - measured electrode
% goodness_mat(el,elm,1) - relative root mean square error
% goodness_mat(el,elm,2) - Spearman's rho
% goodness_mat(el,elm,3) - p-value for Spearman's
%
% goodnessOfFitClassified - el by elm matrix of ones (fitted N1) and zeros(wrong fit)
% Model fits an N1 if:
% Rel RMSE is less than 2.5
% Spearman's rho is larger than .6
% p-value smaller than preset alpha

ccep_parmsN1fit = ccep_setFitParms(t);

[IFG, SMG, STG, els2check] = ccep_anat_locs(s);

goodness_mat = NaN(size(data,1),size(data,1),3);

disp('getting relative error and fit-test correlation data...')
for el = els2check
    if numel(ccep(el).epochs) == 0 % skip this electrode, if it wasn't stimulated.
        disp(['el ' int2str(el) ' not stimulated'])
    else
        for elm = els2check
            epochCheck = sum(squeeze(isnan(data(elm,:,ccep(el).epochs)))); % Check if there are too many bad epochs
            if sum((epochCheck))<length(ccep(el).epochs)/2 % If less than half the epochs are bad, then run model testing
                [errorN1_rel_rms, spearmanMat] = ccep_test_modelChan(data, n1_data, ccep, el, elm, t, ccep_parmsN1fit);
                goodness_mat(el,elm,:) = cat(2,errorN1_rel_rms, spearmanMat);
            else
                continue
            end
        end
    end
end
disp('done')
clear el elm epochCheck rrmse smat

%% In case we've already got the data, and want to only run a part of the script.
% N1 and GoF data should already be loaded.

% switch s
%     case '99'
%         goodnessClassified = goodnessClassified_99;
%         goodness_mat = goodness_mat_99;
%         n1_data = n1_99;
%     case '31'
%         goodnessClassified = goodnessClassified_31;
%         goodness_mat = goodness_mat_31;
%         n1_data = n1_31;
%     case '78'
%         goodnessClassified = goodnessClassified_78;
%         goodness_mat = goodness_mat_78;
%         n1_data = n1_78;
%     case '81'
%         goodnessClassified = goodnessClassified_81;
%         goodness_mat = goodness_mat_81;
%         n1_data = n1_81;
%     case '88'
%         goodnessClassified = goodnessClassified_88;
%         goodness_mat = goodness_mat_88;
%         n1_data = n1_88;
%     case '50'
%         goodnessClassified = goodnessClassified_50;
%         goodness_mat = goodness_mat_50;
%         n1_data = n1_50;
% end

%% Do the classification
%  
%%%%%%%%%%%%%%% We need this data %%%%%%%%%%%%%%%%%%%
s_info = subj_info_ccep(s); % subject info
eval(['ccep = subj_info_SPES_' s]); % CCEP structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set alpha (corrected)
alphaVal = 0.05;
disp('checking goodnes of fit...')
goodnessOfFitClassified = zeros(size(n1_data,2),size(n1_data,2));
for el = els2check
    if numel(ccep(el).epochs) == 0 % skip this electrode, if it wasn't stimulated.
        disp(['el ' int2str(el) ' not stimulated'])
    else
        for elm = els2check
            
                if goodness_mat(el,elm,1) < 2.5 && goodness_mat(el,elm,2) > 0.6 ...
                    && goodness_mat(el,elm,3) < alphaVal                
                goodnessOfFitClassified(el,elm) = 1;
                else
                    continue
                end
        end
    end
end
disp('done')

%% Plot and save goodness of fit matrix

firstEl = min(els2check); 
lastEl = max(els2check); 
figure('Color', 'k')
imagesc(goodnessOfFitClassified)
xlabel('measured')
ylabel('stimulated')
axis([firstEl lastEl firstEl lastEl])
colormap(pinco)
title(['subject ' s ' classified goodness of fit'])
set(gcf, 'InvertHardCopy', 'off');
print(gcf, '-painters', '-dpng', '-r0', strcat(['figures/goodness_fit_mat_s' s]))

clear firstEl lastEl 
% close all
%% Save data

switch s
    case '99'
        goodness_mat_99 = goodness_mat;
        goodnessClassified_99 = goodnessOfFitClassified;
        save('m-files/detectN1/data/goodness99', 'goodness_mat_99', 'goodnessClassified_99');
        clear goodness_mat_99 goodnessClassified_99
        
    case '31'
        goodness_mat_31 = goodness_mat;
        goodnessClassified_31 = goodnessOfFitClassified;
        save('m-files/detectN1/data/goodness31', 'goodness_mat_31', 'goodnessClassified_31');
        clear goodness_mat_31 goodnessClassified_31
        
    case '78'
        goodness_mat_78 = goodness_mat;
        goodnessClassified_78 = goodnessOfFitClassified;
        save('m-files/detectN1/data/goodness78', 'goodness_mat_78', 'goodnessClassified_78');
        clear goodness_mat_78 goodnessClassified_78
           
    case '81'
        goodness_mat_81 = goodness_mat;
        goodnessClassified_81 = goodnessOfFitClassified;
        save('m-files/detectN1/data/goodness81', 'goodness_mat_81', 'goodnessClassified_81');
        clear goodness_mat_81  goodnessClassified_81
                   
    case '88'
        goodness_mat_88 = goodness_mat;
        goodnessClassified_88 = goodnessOfFitClassified;
        save('m-files/detectN1/data/goodness88', 'goodness_mat_88', 'goodnessClassified_88');
        clear goodness_mat_88 goodnessClassified_88
                   
    case '50'
        goodness_mat_50 = goodness_mat;
        goodnessClassified_50 = goodnessOfFitClassified;
        save('m-files/detectN1/data/goodness50', 'goodness_mat_50', 'goodnessClassified_50');
        clear goodness_mat_50 goodnessClassified_50
end
disp('goodness of fit data saved')
