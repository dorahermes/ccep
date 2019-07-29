function [SPESconfig] = ccep_detectERs(stimulus, keepelek, data, fs, numstim)
%
% this function detects early responses after Single Pulse Electrical
% Stimulation
% Inputs:
%
% patconfig contains the electrodes that are in this ECoG file
% spesconfig contains stimulus, which is a struct [the number of stimulus
% pairs x1], with the samplenumbers [10x1] in spesconfig.stimulus.startsamp.
% Detected ERs are saved in spesconfig.stimulus.detERs.
%
% Outputs:
%
% This code adapts early response detection from Dorien van Blooijs to work
% on iEEG BIDS data
%
% original author: Dorien van Blooijs, January 2018
% edited by: Jaap van der Aar, Dora Hermes, Dorien van Blooijs, Giulio
% Castegnaro, UMC Utrecht, 2019


% for each stimulus % I think this is each stimulation pair?
for i=1:size(stimulus,2)
    
    for j=1:size(keepelek,2) % electrodes
        for k=1:numstim % % I think this is the number of times each electrode pair is stimulated: 5 or 10?
            SPES{i,j}(k,:) = data(j,stimulus(i).startsamp(k)-2*fs: stimulus(i).startsamp(k)+2*fs);
        end
        % putting the average across stimuli back in the structure:
        SPES_av(i).smep(j,:) = mean(SPES{i,j});
        % i is a stimulation pair
        % j is a measured electrode
    end
end


% pre-allocation: variables determined in my thesis
thresh = 2.5;
minSD = 50; 
sel = 20;


SDfactor=[];
SD = [];
% pre-allocation
extrasamps = 20; % after the stimulus artefact, no data can be read in 20 samples

for j = 1:size(SPES_av,2)
    smep = SPES_av(j).smep;
    stimchannr = stimulus(j).stimnum;

ERs = [];

for i=1:size(smep,1) % for each electrode
    
    % determine median en sd
    smepmediantotal = median(smep(i,:));
    signal_new = smep(i,:) - smepmediantotal;
    smeprmsbefore = std((signal_new(1,1:round(fs*1.9))));
    
    % when sd < minimal sd, than sd = minimal sd
    if smeprmsbefore < minSD
        smeprmsbefore = minSD;
    end
    
    % when the electrode is stimulated
    if ismember(i,stimchannr)
        sample = 0;
        ampl = 0;
        
        % in other electrode
    elseif ~ismember(i,stimchannr)
        % find positive and negative peak
        [samppos, amplpos]  = peakfinder(signal_new(1,fs*2+extrasamps:round(fs*2.1)),sel,[],1);
        [sampneg, amplneg] = peakfinder(signal_new(1,fs*2+extrasamps:round(fs*2.1)),sel,[],-1);
        
        % excluding the first and last sample
        amplpos(samppos==1) = [];
        samppos(samppos==1) =[];
        amplneg(sampneg==1) = [];
        sampneg(sampneg==1) = [];
        amplpos(samppos >= round(fs*0.1)-extrasamps) = [];
        samppos(samppos >= round(fs*0.1)-extrasamps) = [];
        amplneg(sampneg >= round(fs*0.1)-extrasamps) = [];
        sampneg(sampneg >= round(fs*0.1)-extrasamps) = [];
        
        sampleneg =[];
        amplineg =[];
        if ~isempty(sampneg)
            maxamplneg = find(abs(amplneg) == max(abs(amplneg)));
            sampleneg = sampneg(maxamplneg(1))+fs*2+extrasamps;
            amplineg = amplneg(maxamplneg(1));
        end
        
        samplepos =[];
        amplipos =[];
        if ~isempty(samppos)
            maxamplpos = find(abs(amplpos) == max(abs(amplpos)));
            samplepos = samppos(maxamplpos(1)) +fs*2+extrasamps; % the highest and only sample positive
            amplipos = amplpos(maxamplpos(1));           % the highest and only amplitude positive
        end
        
        % determine the amplitude of the largest peak (positive or negative)
        if ~isempty(samplepos) && isempty(sampleneg) % only samppos
            sample = samplepos;
            ampl = amplipos;
        elseif isempty(samplepos) && ~isempty(sampleneg) % only sampneg
            sample = sampleneg;
            ampl = amplineg;
        elseif ~isempty(samplepos) && ~isempty(sampleneg) % both
            if abs(amplipos) > abs(amplineg) % pos>neg
                sample = samplepos;
                ampl = amplipos;
            elseif abs(amplineg) > abs(amplipos) % neg>pos
                sample = sampleneg;
                ampl = amplineg;
            elseif abs(amplineg) == abs(amplipos) % pos==neg
                if sampleneg < samplepos % pick location of peak
                    sample = sampleneg;
                    ampl = amplineg;
                elseif sampleneg > samplepos
                    sample = samplepos;
                    ampl = amplipos;
                end
            end
        elseif isempty(samplepos) && isempty(sampleneg)
            sample = 0;
            ampl = 0;
        end
    end
    
    % when peak amplitude is saturated, it is deleted
    if abs(ampl) > 3000
        ampl = 0;
    end
    
    % is a peak an ER or not?
    if abs(ampl) > thresh* abs(smeprmsbefore)
        ERs = [ERs i];
    end
    
    SDfactor(i) = abs(ampl)/abs(smeprmsbefore);
    SD(i) = smeprmsbefore;
    Amplitude(i) = ampl;
end

stimulus(j).detERs = ERs; % ERs for all stimulation pairs

end

SPESconfig.stimulus = stimulus;

end


