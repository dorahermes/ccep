%% Get connectivity data for subject
subjects = {'99','31','78','81','88'};
s_nr = 2;

s = subjects{s_nr};
eval(['ccep = subj_info_SPES_' s]); % make CCEP structure
s_info = subj_info_ccep(s); % get subject info
load(s_info.cortex) % get cortex
load(s_info.els) 

switch s_nr
    case 1
        adjAmp = adjAmp_99;
        recData = recData_99;
        n1_data = n1_99;
        projData = projData_99;
    case 2
        adjAmp = adjAmp_31;
        recData = recData_31;
        n1_data = n1_31;
        projData = projData_31;
    case 3
        adjAmp = adjAmp_78;
        recData = recData_78;
        n1_data = n1_78;
        projData = projData_78;
    case 4
        adjAmp = adjAmp_81;
        recData = recData_81;
        n1_data = n1_81;
        projData = projData_81;
    case 5
        adjAmp = adjAmp_88;
        recData = recData_88;
        n1_data = n1_88;
        projData = projData_88;
end
