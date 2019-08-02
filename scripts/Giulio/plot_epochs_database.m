
%Warning: it is best to run workflow and
%database_ccep_visual_conn_matrix_N1.m before this script 

%% ADAPT THIS CODE TO THE DATABASE TO LOOK AT THE AVERAGED EPOCHS 
%% plot averaged epochs

tt = [1:database(subj).metadata(runs).epoch_length*database(subj).metadata(runs).data_hdr.Fs]/database(subj).metadata(runs).data_hdr.Fs - database(subj).metadata(runs).epoch_prestim_length; 


stim_pair_nr = 1; 
ccep_elec = 3; 
figure()
plot(tt(tt>-.5 & tt<0.2),zeros(size(tt(tt>-.5 & tt<0.2))),'Color',[.5 .5 .5])

%plot(tt,squeeze(cc_epoch_sorted_avg(ccep_elec,stim_pair_nr,:)))
plot(tt(tt>-.5 & tt<0.2),squeeze(database(subj).metadata(runs).epoched_data_avg(ccep_elec,stim_pair_nr,(tt>-.5 & tt<0.2))))
xlabel('time(s)')
% set(gca,'Ydir','reverse')
ylabel('amplitude(uV)')
title(['elec ' database(subj).metadata(runs).electrodes.name(ccep_elec) ' for stimulation of ' database(subj).metadata(runs).electrodes.name(database(subj).metadata(runs).stimulated_pairs(stim_pair_nr,1)) ' and ' database(subj).metadata(runs).electrodes.name(database(subj).metadata(runs).stimulated_pairs(stim_pair_nr,2))])


%% THIS IS TO PLOT ALL THE SIGNALS FOR ONE STIMULATION 
%% plot to test
extrasamps = 5; % set to 5 now, to read properties of N1 onset
baseline_tt = tt>-2 & tt<-.1;
% makes it possible to plot multiple CCEPs to check
for ii = 3 %decide the electrode to look at 
    for jj = 1:56 % decide the stimulation to look at 
        
        % recalculate new_signal for this ii and jj to plot
        signal_median = median(database(subj).metadata(runs).epoched_data_avg(ii,jj,baseline_tt),3);
        new_signal = squeeze(database(subj).metadata(runs).epoched_data_avg(ii,jj,:)) - signal_median; 
    
        plot((find(tt>0,1)+extrasamps:find(tt>0.1,1)),new_signal(find(tt>0,1)+extrasamps:find(tt>0.1,1)))
        hold on

        plot(database(subj).metadata(runs).n1_peak_sample(ii,jj),database(subj).metadata(runs).n1_peak_amplitude(ii,jj),'r*')
    end
end 



%% plot signals in all electrodes with a label 
figure()

stim_pair_nr = 20; % decide stimulated pair
%database(subj).metadata(runs).labeled_stimsets run this line to visualize
%which pair is stimulated rapidly 

suptitle(['Electrode by electrode analysis of the stimulation of' database(subj).metadata(runs).electrodes.name(database(subj).metadata(runs).stimulated_pairs(stim_pair_nr,1)),...
    ' to ' database(subj).metadata(runs).electrodes.name(database(subj).metadata(runs).stimulated_pairs(stim_pair_nr,2))])

% plot all measured channels
nr_column = 10;
nr_rows = ceil(length(database(subj).metadata(runs).channels.name)/nr_column);

% makes it possible to plot multiple CCEPs to check (plot only signals for
% which there is a new_label assigned to it)

for ii = 1:length(database(subj).metadata(runs).channels.name)
    %if ~isnan(database(subj).metadata(runs).electrodes.Wang_label(ii)) 
    if isequal(database(subj).metadata(runs).channels.type(ii),{'ECOG'}) && ismember(database(subj).metadata(runs).new_label(ii),[1:25]) 
        
        subplot(nr_rows,nr_column,ii);
        hold on
        title(['Area' num2cell(database(subj).metadata(runs).new_label(ii))])
        
        % recalculate new_signal for this ii and jj to plot
        signal_median = median(database(subj).metadata(runs).epoched_data_avg(ii,stim_pair_nr,baseline_tt),3);
        new_signal = squeeze(database(subj).metadata(runs).epoched_data_avg(ii,stim_pair_nr,:)) - signal_median; 
        plot(tt(find(tt>0,1)+extrasamps:find(tt>0.1,1)),new_signal(find(tt>0,1)+extrasamps:find(tt>0.1,1)))

        hold on
        if ~isnan(database(subj).metadata(runs).n1_peak_sample(ii,stim_pair_nr)) % if a N1 is detected
            plot(tt(database(subj).metadata(runs).n1_peak_sample(ii,stim_pair_nr)),database(subj).metadata(runs).n1_peak_amplitude(ii,stim_pair_nr),'r*')
        end 
    end 
end

