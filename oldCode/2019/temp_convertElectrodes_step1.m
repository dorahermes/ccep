% load the electrodes for a subject
load('/Fridge/bci/data/14-420_adults/franeker/analysed/ALICE/results/projected_electrodes_coord/franeker_projectedElectrodes_FreeSurfer_3dclust.mat')


% save it as temp
save('/Fridge/users/jaap/ccep/dataBIDS/sourcedata/temp/temp.mat','elecmatrix')
