% [DATAOUT]=readtrc(SETTIINGS)  -   Reads Micromed System 98 *.trc EEG File
%                                   for import into EEGLAB. Tghi file reads
%                                   Micromed System Plus EEG *.trc files with
%                                   header of type 4
%
% USAGE:
%   >> [DATAOUT]=readtrc(SETTINGS)
%
%
% INPUT
%   SETTINGS is a struct holding the parameters for reading the .TRC file
%   SETTINGS has the following fields:
%
%       SETTINGS.filename :                 Name of file to be imported
%       SETTINGS.loadevents.state :         'yes' for loading event triggers
%                                           'no' for not
%                                           default=yes
%       SETTINGS.loadevents.type :          'analog' for event triggers inserted
%                                           on 'MKR channel
%                                           'digital' for triggers inserted on
%                                           EEG channels
%                                           default=analog
%       SETTINGS.loadevents.dig_ch1:        number of name of digital channel 1
%       SETTINGS.loadevents.dig_ch1_label:  label to give events on digital channel 1
%       SETTINGS.loadevents.dig_ch2:        number of name of digital channel 2
%       SETTINGS.loadevents.dig_ch2_label:  label to give events on digital channel 2
%       SETTINGS.chan_adjust_status:        1 for adjusting amp of channels 0 for not
%                                           default=0
%       SETTINGS.chan_exclude_status        1 for excluding some channels 0 for not
%                                           default=no
%       SETTINGS.chan_adjust                channels to adjust
%       SETTINGS.chan_exclude               channels to exclude
%       SETTINGS.verbose                    output extra fields: title,
%                                           laboratory, birth_month, birth_day,
%                                           acquisition_unit,
%                                           acquisition_unit_description,
%                                           filetype,filetype_code,filetype_description,
%                                           Code_Area,Electrode_Area,
%                                           electrode,originalelectrode
%                                           default=no
%       SETTINGS.loaddata.state             'yes' for loading data
%                                           'no' for only header info
%                                           default=yes
%       SETTINGS.loaddata.type              'uV' for EEGLAB
%                                           'raw'for comparison with BCI2000
%       SETTINGS.movenonEEGchannels         moves markerchannels and oxymeter channels
%                                           to last channels (e.g. with multiple LTM
%                                           headboxes)
%                                           default='yes'
%
%   Alternat method: enter 11 inputr arguments each corresponding to one of the above
%   fields.  If only one arg is used,  it must be the struct above.  If more, there
%   must be 11 inputs in the order above i.e. OUT=readtrc(filename,eventstate....etc).
%
% OUTPUT
%   DATAOUT with same fields as EEGLAB EEG file structure. (see eeg_checkset.m)
%   The following fields are use (with verbose='no'):
%       DATAOUT.data
%       DATAOUT.filename
%       DATAOUT.filepath
%       DATAOUT.srate
%       DATAOUT.setname
%       DATAOUT.pnts
%       DATAOUT.nbchan
%       DATAOUT.trials
%       DATAOUT.xmin
%       DATAOUT.ref
%
% Copyright (c) 2004 Rami Niazy, FMRIB Centre, University of Oxford.
% Additions (c) 2013 Erik Aarnoutse, UMC Utrecht.
%
%   see also pop_readtrc.m    eegplugin_trcimport.m

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Rami K. Niazy, FMRIB Centre, University of Oxford
%                    niazy@ieee.org
%
% This program can only be provided from Micromed s.r.l.  You may not
% give away or edit this program.

% Sep 9, 2013
% Added markerchannels move to end of channel list

% Dec 3, 2012
% Added verbose output

% Sep 29, 2004
% Fixed Trigger offset

% Sep 22, 2004
% Fixed Trigger offset

% Sep 21, 2004
% Allow usage of ':' in sepcifying
% channels to exclude or adjust

% Aug 12, 2004
% Added command line usage
% fixed com out for EEGLAB

% Aug 10, 2004
% Adjusted Licence
% Added Chan exclude
% Fixed scaling

% Aug 4, 2004
% Fixed scaling of data to uV
% Fixed Analog Trig size detection

% Date: Julyl 30, 2004
% Initial Setup

function [TRC]=jun_readtrc(varargin)

if nargin<1
    error('Not enough input arguments. Please see help file for usage');
elseif nargin==1
    if isstruct(varargin{1})
        PARAM=varargin{1};
    else
        error('When using single input argument, it must be a Structure.  Please see help file.');
    end
elseif nargin < 11
    error('Not enough input arguments. Please see help file for usage');
elseif nargin > 11
    error('Too many input arguments.  Please see help file for usage');
else
    PARAM.filename=varargin{1};
    PARAM.loadevents.state=varargin{2};
    PARAM.loadevents.type=varargin{3};
    PARAM.loadevents.dig_ch1=varargin{4};
    PARAM.loadevents.dig_ch1_label=varargin{5};
    PARAM.loadevents.dig_ch2=varargin{6};
    PARAM.loadevents.dig_ch2_label=varargin{7};
    PARAM.chan_adjust_status=varargin{8};
    PARAM.chan_exclude_status=varargin{9};
    PARAM.chan_adjust=varargin{10};
    PARAM.chan_exclude=varargin{11};
end

% defaults added by E.J. Aarnoutse
if ~isfield(PARAM,'verbose')
    PARAM.verbose='no';
end
if ~isfield(PARAM,'loadevents')
    PARAM.loadevents.state='yes';
    PARAM.loadevents.type='analog';
end
if ~isfield(PARAM,'chan_adjust_status')
    PARAM.chan_adjust_status=0;
end
if ~isfield(PARAM,'chan_exclude_status')
    PARAM.chan_exclude_status=0;
end
if ~isfield(PARAM,'loaddata')
    PARAM.loaddata.state='yes';
end
if ~isfield(PARAM.loaddata,'type')
    PARAM.loaddata.type='uV';
end
if ~isfield(PARAM,'movenonEEGchannels')
    PARAM.movenonEEGchannels='yes';
end


% ---------------- Opening File------------------
trcfile=PARAM.filename;
fid=fopen(trcfile,'r','l','US-ASCII');
if fid==-1
    error('Can''t open *.trc file')
end
%TRC=eeg_emptyset;
fprintf('Reading %s...\n',trcfile);

%------------------reading patient & recording info----------
trctitle=char(fread(fid,32,'char'))';         % offset 0
laboratory=char(fread(fid,32,'char'))';       % offset 32
%fseek(fid,64,-1);
surname=char(fread(fid,22,'char'))';          % offset 64
name=char(fread(fid,20,'char'))';             % offset 76
birth_month=fread(fid,1,'uint8');             % offset 96
birth_day=fread(fid,1,'uint8');               % offset 97
birth_year=1900 + fread(fid,1,'uint8');       % offset 98

fseek(fid,128,-1);
day=fread(fid,1,'char');                      % offset 128
if length(num2str(day))<2
    day=['0' num2str(day)];
else
    day=num2str(day);
end
monthnr=fread(fid,1,'char');                  % offset 129
switch monthnr
    case 1
        month='JAN';
    case 2
        month='FEB';
    case 3
        month='MAR';
    case 4
        month='APR';
    case 5
        month='MAY';
    case 6
        month='JUN';
    case 7
        month='JUL';
    case 8
        month='AUG';
    case 9
        month='SEP';
    case 10
        month='OCT';
    case 11
        month='NOV';
    case 12
        month='DEC';
    otherwise
        month = '';
end
year=num2str(fread(fid,1,'char')+1900);       % offset 130
hour=fread(fid,1,'char');                     % offset 131
minutes=fread(fid,1,'char');                  % offset 132
seconds=fread(fid,1,'char');                  % offset 133

%------------------ Reading Header Info ---------
acquisition_unit=fread(fid,1,'uint16');        % offset 134
switch acquisition_unit
    case 0
        acquisition_unit_description='BQ124 - 24 channels headbox, Internal Interface';
    case 1
        acquisition_unit_description='MS40 - Holter recorder';
    case 2
        acquisition_unit_description='BQ132S - 32 channels headbox, Internal Interface';
    case 6
        acquisition_unit_description='BQ124 - 24 channels headbox, BQ CARD Interface';
    case 7
        acquisition_unit_description='SAM32 - 32 channels headbox, BQ CARD Interface';
    case 8
        acquisition_unit_description='SAM25 - 25 channels headbox, BQ CARD Interface';
    case 9
        acquisition_unit_description='BQ132S R - 32 channels reverse headbox, Internal Interface';
    case 10
        acquisition_unit_description='SAM32 R - 32 channels reverse headbox, BQ CARD Interface';
    case 11
        acquisition_unit_description='SAM25 R - 25 channels reverse headbox, BQ CARD Interface';
    case 12
        acquisition_unit_description='SAM32 - 32 channels headbox, Internal Interface';
    case 13
        acquisition_unit_description='SAM25 - 25 channels headbox, Internal Interface';
    case 14
        acquisition_unit_description='SAM32 R - 32 channels reverse headbox, Internal Interface';
    case 15
        acquisition_unit_description='SAM25 R - 25 channels reverse headbox, Internal Interface';
    case 16
        acquisition_unit_description='SD - 32 channels headbox with jackbox, SD CARD Interface ? PCI Internal Interface';
    case 17
        acquisition_unit_description='SD128 - 128 channels headbox, SD CARD Interface ? PCI Internal Interface';
    case 18
        acquisition_unit_description='SD96 - 96 channels headbox, SD CARD Interface ? PCI Internal Interface';
    case 19
        acquisition_unit_description='SD64 - 64 channels headbox, SD CARD Interface ? PCI Internal Interface';
    case 20
        acquisition_unit_description='SD128c - 128 channels headbox with jackbox, SD CARD Interface ? PCI Internal Interface';
    case 21
        acquisition_unit_description='SD64c - 64 channels headbox with jackbox, SD CARD Interface ? PCI Internal Interface';
    case 22
        acquisition_unit_description='BQ132S - 32 channels headbox, PCI Internal Interface';
    case 23
        acquisition_unit_description='BQ132S R - 32 channels reverse headbox, PCI Internal Interface';
    case 49
        acquisition_unit_description='LTM express 256 channels headbox,  ? PCI Internal Interface';
    otherwise
        acquisition_unit_description='unknown acquisition_unit';
end

filetype=fread(fid,1,'uint16');                % offset 136
switch filetype
    
    case 40
        filetype_code='C128';
        filetype_description='C.R., 128 EEG (headbox SD128 only)';
    case 42
        filetype_code='C84P';
        filetype_description='C.R., 84 EEG, 44 poly (headbox SD128 only)';
    case 44
        filetype_code='C84';
        filetype_description='C.R., 84 EEG, 4 reference signals (named MKR,MKRB,MKRC,MKRD) (headbox SD128 only)';
    case 46
        filetype_code='C96';
        filetype_description='C.R., 96 EEG (headbox SD128 ? SD96 ? BQ123S(r))';
    case 48
        filetype_code='C63P';
        filetype_description='C.R., 63 EEG, 33 poly';
    case 50
        filetype_code='C63';
        filetype_description='C.R., 63 EEG, 3 reference signals (named MKR,MKRB,MKRC)';
    case 52
        filetype_code='C64';
        filetype_description='C.R., 64 EEG';
    case 54
        filetype_code='C42P';
        filetype_description='C.R., 42 EEG, 22 poly';
    case 56
        filetype_code='C42';
        filetype_description='C.R., 42 EEG, 2 reference signals (named MKR,MKRB)';
    case 58
        filetype_code='C32';
        filetype_description='C.R., 32 EEG';
    case 60
        filetype_code='C21P';
        filetype_description='C.R., 21 EEG, 11 poly';
    case 62
        filetype_code='C21';
        filetype_description='C.R., 21 EEG, 1 reference signal (named MKR)';
    case 64
        filetype_code='C19P';
        filetype_description='C.R., 19 EEG, variable poly';
    case 66
        filetype_code='C19';
        filetype_description='C.R., 19 EEG, 1 reference signal (named MKR)';
    case 68
        filetype_code='C12';
        filetype_description='C.R., 12 EEG';
    case 70
        filetype_code='C8P';
        filetype_description='C.R., 8 EEG, variable poly';
    case 72
        filetype_code='C8';
        filetype_description='C.R., 8 EEG';
    case 74
        filetype_code='CFRE';
        filetype_description='C.R., variable EEG, variable poly';
    case 76
        filetype_code='C25P';
        filetype_description='C.R., 25 EEG (21 standard, 4 poly transformed to EEG channels), 7 poly ? headbox BQ132S(r) only';
    case 78
        filetype_code='C27P';
        filetype_description='C.R., 27 EEG (21 standard, 6 poly transformed to EEG channels), 5 poly ? headbox BQ132S(r) only';
    case 80
        filetype_code='C24P';
        filetype_description='C.R., 24 EEG (21 standard, 3 poly transformed to EEG channels), 8 poly ? headbox SAM32(r) only';
    case 82
        filetype_code='C25P';
        filetype_description='C.R., 25 EEG (21 standard, 4 poly transformed to EEG channels), 7 poly ? headbox SD with headbox JB 21P';
    case 84
        filetype_code='C27P';
        filetype_description='C.R., 27 EEG (21 standard, 6 poly transformed to EEG channels), 5 poly ? headbox SD with headbox JB 21P';
    case 86
        filetype_code='C31P';
        filetype_description='C.R., 27 EEG (21 standard, 10 poly transformed to EEG channels), 1 poly ? headbox SD with headbox JB 21P';
    case 100
        filetype_code='C26P';
        filetype_description='C.R., 26 EEG, 6 poly (headbox SD, SD64c, SD128c with headbox JB Mini)';
    case 101
        filetype_code='C16P';
        filetype_description='C.R., 16 EEG, 16 poly (headbox SD with headbox JB M12)';
    case 102
        filetype_code='C12P';
        filetype_description='C.R., 12 EEG, 20 poly (headbox SD with headbox JB M12)';
    case 103
        filetype_code='32P';
        filetype_description='32 poly (headbox SD, SD64c, SD128c with headbox JB Bip)';
    case 120
        filetype_code='C48P';
        filetype_description='C.R., 48 EEG, 16 poly (headbox SD64)';
    case 121
        filetype_code='C56P';
        filetype_description='C.R., 56 EEG, 8 poly (headbox SD64)';
    case 122
        filetype_code='C24P';
        filetype_description='C.R., 24 EEG, 8 poly (headbox SD64)';
    case 140
        filetype_code='C52P';
        filetype_description='C.R., 52 EEG, 12 poly (headbox SD64c, SD128c with 2 headboxes JB Mini)';
    case 141
        filetype_code='64P';
        filetype_description='64 poly (headbox SD64c, SD128c with 2 headboxes JB Bip)';
    case 160
        filetype_code='C88P';
        filetype_description='C.R., 88 EEG, 8 poly (headbox SD96)';
    case 161
        filetype_code='C80P';
        filetype_description='C.R., 80 EEG, 16 poly (headbox SD96)';
    case 162
        filetype_code='C72P';
        filetype_description='C.R., 72 EEG, 24 poly (headbox SD96)';
    case 180
        filetype_code='C120P';
        filetype_description='C.R., 120 EEG, 8 poly (headbox SD128)';
    case 181
        filetype_code='C112P';
        filetype_description='C.R., 112 EEG, 16 poly (headbox SD128)';
    case 182
        filetype_code='C104P';
        filetype_description='C.R., 104 EEG, 24 poly (headbox SD128)';
    case 183
        filetype_code='C96P';
        filetype_description='C.R., 96 EEG, 32 poly (headbox SD128)';
    case 200
        filetype_code='C122P';
        filetype_description='C.R., 122 EEG, 6 poly (headbox SD128c with 4 headboxes JB Mini)';
    case 201
        filetype_code='C116P';
        filetype_description='C.R., 116 EEG, 12 poly (headbox SD128c with 4 headboxes JB Mini)';
    case 202
        filetype_code='C110P';
        filetype_description='C.R., 110 EEG, 18 poly (headbox SD128c with 4 headboxes JB Mini)';
    case 203
        filetype_code='C104P';
        filetype_description='C.R., 104 EEG, 24 poly (headbox SD128c with 4 headboxes JB Mini)';
    case 204
        filetype_code='128P';
        filetype_description='128 poly (headbox SD128c with 4 headboxes JB Bip)';
    case 205
        filetype_code='96P';
        filetype_description='96 poly (headbox SD128c with 3 headboxes JB Bip)';
    otherwise
        filetype_code='';
        filetype_description='unknown filetype';
end



fseek(fid,175,-1);
Header_Type=fread(fid,1,'char');
if Header_Type ~= 4
    error('*.trc file is not Micromed System98 Header type 4')
end

fseek(fid,138,-1);
Data_Start_Offset=fread(fid,1,'uint32');       % offset 138
Num_Chan=fread(fid,1,'uint16');                % offset 142
Multiplexer=fread(fid,1,'uint16');             % offset 144
Rate_Min=fread(fid,1,'uint16');                % offset 146
Bytes=fread(fid,1,'uint16');                   % offset 148

%Compression
%Montages
%Dvideo_Begin
%MPEG_Delay
%Reserved_1
fseek(fid,176+8,-1);

Code_Area=fread(fid,1,'uint32');               % offset 184
Code_Area_Length=fread(fid,1,'uint32');        % offset 188

fseek(fid,192+8,-1);

Electrode_Area=fread(fid,1,'uint32');          % offset 200
Electrode_Area_Length=fread(fid,1,'uint32');   % offset 204

fseek(fid,400+8,-1);
Trigger_Area=fread(fid,1,'uint32');
Trigger_Area_Length=fread(fid,1,'uint32');

%------------ Read Channel information ------
fseek(fid,1154,-1);
fseek(fid,266738,-1);
fseek(fid,328178,-1);

nrchannels=Num_Chan;
for ch=1:nrchannels
    channelx(ch).name = sprintf('%s',fread(fid,6,'char'));
    channelx(ch).ref = char(fread(fid,6,'char'));
    channelx(ch).unknown = fread(fid,116,'char');
end

fseek(fid,5250,-1);

nrchannels=Num_Chan;
for ch=1:nrchannels
    channel(ch).name = sprintf('%s',fread(fid,6,'char'));
    channel(ch).ref = sprintf('%s',fread(fid,6,'char'));
    channel(ch).unknown = fread(fid,116,'char')';
end

%----------------- Read Trace Data ----------

if strcmp(PARAM.loaddata.state,'yes')
    
    fseek(fid,Data_Start_Offset,-1);
    switch Bytes
        case 1
            mmtrace=fread(fid,'uint8');
        case 2
            mmtrace=fread(fid,'uint16');
        case 4
            mmtrace=fread(fid,'uint32');
            %disp('data32read');
    end
    
    m=length(mmtrace);
    if rem(m,Num_Chan)~=0
        roundata=floor(m/Num_Chan);
        mmtrace=mmtrace(1:roundata*Num_Chan);
        m=length(mmtrace);
    end
    mmtrace=reshape(mmtrace,Num_Chan,m/Num_Chan);
    m=length(mmtrace);
else
    fseek(fid,0,'eof');
    e=ftell(fid);
    m=floor((e-Data_Start_Offset)/(Bytes*Num_Chan));
    
end
%------------------ Reading Code Info -------------
fseek(fid,Code_Area,-1);
code=fread(fid,Num_Chan,'uint16');


for c=1:Num_Chan
    electrode(c).chan_record=code(c);
    fseek(fid,Electrode_Area+code(c)*128,-1);
    electrode(c).status=fread(fid,1,'uint8');
    electrode(c).type_code=fread(fid,1,'uint8');
    
    if (bitget(electrode(c).type_code,1)==0)
        electrode(c).reference = 'G2';
    else
        electrode(c).reference = 'bipolar';
    end
    electrode(c).marker = bitget(electrode(c).type_code,2);
    electrode(c).oxym = bitget(electrode(c).type_code,3);
    electrode(c).t_16DC = bitget(electrode(c).type_code,4);
    electrode(c).bip2eeg = bitget(electrode(c).type_code,5);
    
    electrode(c).label=char(fread(fid,6,'char'))';
    if c <10
        electrode(c).positive_input=[num2str(c),' -',electrode(c).label];
        electrode(c).negative_input=[num2str(c),' -',char(fread(fid,6,'char'))'];
    else
        electrode(c).positive_input=[num2str(c),'-',electrode(c).label];
        electrode(c).negative_input=[num2str(c),'-',char(fread(fid,6,'char'))'];
    end
    
    electrode(c).logical_min=fread(fid,1,'int32');
    
    electrode(c).logical_max=fread(fid,1,'int32');
    
    electrode(c).logical_ground=fread(fid,1,'int32');
    
    electrode(c).physical_min=fread(fid,1,'int32');
    
    electrode(c).physical_max=fread(fid,1,'int32');
    
    electrode(c).measurement_unit=fread(fid,1,'int16');
    
    switch electrode(c).measurement_unit
        case -1
            electrode(c).measurement_unit=1e-9;
        case 0
            electrode(c).measurement_unit=1e-6;
        case 1
            electrode(c).measurement_unit=1e-3;
        case 2
            electrode(c).measurement_unit=1;
        case 100
            electrode(c).measurement_unit='percent';
        case 101
            electrode(c).measurement_unit='bpm';
        case 102
            electrode(c).measurement_unit='Adim';
        otherwise
            warning('Unknown measurement unit. uV assumed.');
            electrode(c).measurement_unit=10e-6;
    end
    
    electrode(c).prefiltering_hipass_limit=double(fread(fid,1,'int16'))/1000;
    electrode(c).prefiltering_hipass_type=fread(fid,1,'int16');
    electrode(c).prefiltering_lowPass_limit=fread(fid,1,'int16');
    electrode(c).prefiltering_lowPass_type=fread(fid,1,'int16');
    electrode(c).rate_coef=fread(fid,1,'uint16');
    %Sampling Rate Coefficient
    %1 ? min. Sampling Rate
    %2 ? 2?min. Sampling Rate 4 ? 4?min. Sampling Rate
    electrode(c).position=fread(fid,1,'int16');
    electrode(c).latitudine=fread(fid,1,'float');
    electrode(c).longitudine=fread(fid,1,'float');
    electrode(c).presentinmap=fread(fid,1,'uint8');
    electrode(c).isinavg=fread(fid,1,'uint8');
    electrode(c).description=fread(fid,32,'*char')';
    electrode(c).x=fread(fid,1,'float');
    electrode(c).y=fread(fid,1,'float');
    electrode(c).z=fread(fid,1,'float');
    electrode(c).coordinate_type=fread(fid,1,'int16');
    
    if ~strcmpi(PARAM.loaddata.type,'raw') && strcmpi(PARAM.loaddata.state,'yes')
        if ischar(electrode(c).measurement_unit)==0
            mmtrace(c,:)=-((mmtrace(c,:)-electrode(c).logical_ground)/...
                (electrode(c).logical_max-electrode(c).logical_min+1))*...
                (electrode(c).physical_max-electrode(c).physical_min)*electrode(c).measurement_unit;
        else
            mmtrace(c,:)=-((mmtrace(c,:)-electrode(c).logical_ground)/...
                (electrode(c).logical_max-electrode(c).logical_min+1))*...
                (electrode(c).physical_max-electrode(c).physical_min);
        end
    end
end

%---------------- Move Non-EEG channels to end ----------

originalelectrode=electrode;
if strcmpi(PARAM.loaddata.state,'yes')
    numbernoneegchannels=0;
    newc=0;
    c=1;
    while c <= Num_Chan
        if electrode(c).marker==1 || electrode(c).oxym==1
            numbernoneegchannels=numbernoneegchannels+1;
            noneegtrace(numbernoneegchannels,:)=mmtrace(c,:);
            noneegtraces(numbernoneegchannels)=c;
        else
            newc=newc+1;
            newelectrode(newc)=electrode(c);
        end
        c=c+1;
    end
    if strcmpi(PARAM.movenonEEGchannels,'yes') && (numbernoneegchannels > 0)
        mmtrace(noneegtraces,:)=[]; %delete noneegchannels from mmtrace
        mmtrace=vertcat(mmtrace,noneegtrace); %add noneegchannels to end
        for c=1:numbernoneegchannels
            newelectrode(Num_Chan-numbernoneegchannels+c)=electrode(noneegtraces(c));
            if electrode(noneegtraces(c)).marker==1
                disp(['Marker channel ' num2str(noneegtraces(c)) ' moved to channel ' num2str(Num_Chan-numbernoneegchannels+c)]);
            else
                disp(['Oxymeter channel ' num2str(noneegtraces(c)) ' moved to channel ' num2str(Num_Chan-numbernoneegchannels+c)]);
            end
        end
        electrode=newelectrode;
    end
end
%---------------- Reading Trigger Data ----------
fseek(fid,Trigger_Area,-1);
for l=1:Trigger_Area_Length/6
    trigger(1,l)=fread(fid,1,'uint32');
    trigger(2,l)=fread(fid,1,'uint16');
end
%------------------------------------------------

%---------------- Reading Trigger Data ----------
fseek(fid,Electrode_Area,-1);
for l=1:Electrode_Area_Length/128
    
end

%--------------- End reading ----------------------

fclose(fid);

%--------------- Preparing output -----------------
%


first_trigger=trigger(1,1);

tl=length(trigger);
NoTrig=0;
for tr=1:tl
    if ((trigger(1,tr) <= m) && (trigger(1,tr) >= first_trigger))
        NoTrig=NoTrig+1;
    end
end

if NoTrig > 0
    trigger=trigger(:,1:NoTrig);
else
    trigger=[];
    first_trigger=[];
end




mean_fs=mean(cat(1,electrode.rate_coef));
switch mean_fs
    case 0
        fs=1*Rate_Min;
    case 1
        fs=1*Rate_Min;
    case 2
        fs=2*Rate_Min;
    case 3
        fs=3*Rate_Min;
    case 4
        fs=4*Rate_Min;
    case 5
        fs=5*Rate_Min;
    otherwise
        warning('Unsupported Sampling Frequency');
end



%------------- Passing Back Arguments-------------
switch  PARAM.loadevents.state
    case 'no';
    case 'yes';
        TRC.event.type=[];
        TRC.event.latency=[];
        TRC.event.urevent=[];
        TRC.event.value=[];
        TRC.urevent.type=[];
        TRC.urevent.latency=[];
        TRC.urevent.value=[];
        switch PARAM.loadevents.type
            case 'analog'
                fprintf(['Extracting ' num2str(size(trigger,2)) ' events...']);
                if ~isempty(trigger)
                    for E=1:size(trigger,2)
                        TRC.event(E).type='MARKER';
                        TRC.event(E).latency=trigger(1,E)+1;
                        TRC.event(E).urevent=E;
                        TRC.event(E).value=trigger(2,E);
                        TRC.urevent(E).type='MARKER';
                        TRC.urevent(E).latency=trigger(1,E)+1;
                        TRC.urevent(E).value=trigger(2,E);
                    end
                else
                    warndlg('No analog triggers to import on ''MRK'' channel','Import .TRC Warning!');
                end
            case 'digital'
                
                if strcmpi(PARAM.loaddata.state,'yes')
                    fprintf('Extracting events...');
                    %------find 1st trigger channel-------
                    dig_ch=[];
                    if str2num(PARAM.loadevents.dig_ch1)
                        dig_ch(1)=str2num(PARAM.loadevents.dig_ch1);
                    else
                        for C=1:Num_Chan
                            if ~isempty(findstr(lower(PARAM.loadevents.dig_ch1),lower(electrode(C).positive_input)))
                                dig_ch(1)=C;
                                break;
                            end
                        end
                    end
                    
                    %------check for 2nd trigger channel-------
                    if ~isempty(PARAM.loadevents.dig_ch2)
                        if str2num(PARAM.loadevents.dig_ch2)
                            dig_ch(2)=str2num(PARAM.loadevents.dig_ch2);
                        else
                            for C=1:Num_Chan
                                if ~isempty(findstr(lower(PARAM.loadevents.dig_ch2),lower(electrode(C).positive_input)))
                                    dig_ch(2)=C;
                                    break;
                                end
                            end
                        end
                    end
                    
                    
                    %--------read trigs-------------------------
                    Trigs1=[];
                    Trigs2=[];
                    
                    trigval=min(mmtrace(dig_ch(1),:));
                    Trigs1=find(mmtrace(dig_ch(1),:)==trigval);
                    mmtrace(dig_ch(1),Trigs1)=(mmtrace(dig_ch(1),(Trigs1+1))+mmtrace(dig_ch(1),(Trigs1-1)))/2;
                    if length(dig_ch)>1
                        Trigs2=find(mmtrace(dig_ch(2),:)==trigval);
                        mmtrace(dig_ch(2),Trigs2)=(mmtrace(dig_ch(2),(Trigs2+1))+mmtrace(dig_ch(2),(Trigs2-1)))/2;
                    end
                    
                    
                    if length(dig_ch)<2
                        
                        for E=1:length(Trigs1)
                            TRC.event(E).type=PARAM.loadevents.dig_ch1_label;
                            TRC.event(E).latency=Trigs1(E);
                            TRC.event(E).urevent=E;
                            TRC.urevent(E).type=PARAM.loadevents.dig_ch1_label;
                            TRC.urevent(E).latency=Trigs1(E);
                        end
                        
                    else
                        if ~isempty(Trigs2)
                            E=1;
                            E1=1;
                            E1end=0;
                            E2=1;
                            E2end=0;;
                            while ((E1end==0) | (E2end==0))
                                
                                if (E1end==0 & E2end==0)
                                    
                                    if Trigs1(E1) < Trigs2(E2)
                                        TRC.event(E).type=PARAM.loadevents.dig_ch1_label;
                                        TRC.event(E).latency=Trigs1(E1);
                                        TRC.event(E).urevent=E;
                                        TRC.urevent(E).type=PARAM.loadevents.dig_ch1_label;
                                        TRC.urevent(E).latency=Trigs1(E1);
                                        if E1==length(Trigs1)
                                            E1end=1;
                                        end
                                        E=E+1; E1=E1+1;
                                    elseif Trigs1(E1) > Trigs2(E2)
                                        TRC.event(E).type=PARAM.loadevents.dig_ch2_label;
                                        TRC.event(E).latency=Trigs2(E2);
                                        TRC.event(E).urevent=E;
                                        TRC.urevent(E).type=PARAM.loadevents.dig_ch2_label;
                                        TRC.urevent(E).latency=Trigs2(E2);
                                        if E2==length(Trigs2)
                                            E2end=1;
                                        end
                                        E=E+1; E2=E2+1;
                                    end
                                    
                                elseif (E1end==1 & E2end==0)
                                    
                                    if E2==length(Trigs2)
                                        E2end=1;
                                    end
                                    
                                    TRC.event(E).type=PARAM.loadevents.dig_ch2_label;
                                    TRC.event(E).latency=Trigs2(E2);
                                    TRC.event(E).urevent=E;
                                    TRC.urevent(E).type=PARAM.loadevents.dig_ch2_label;
                                    TRC.urevent(E).latency=Trigs2(E2);
                                    E=E+1; E2=E2+1;
                                    
                                elseif (E1end==0 & E2end==1)
                                    
                                    if E1==length(Trigs1)
                                        E1end=1;
                                    end
                                    
                                    TRC.event(E).type=PARAM.loadevents.dig_ch1_label;
                                    TRC.event(E).latency=Trigs1(E1);
                                    TRC.event(E).urevent=E;
                                    TRC.urevent(E).type=PARAM.loadevents.dig_ch1_label;
                                    TRC.urevent(E).latency=Trigs1(E1);
                                    E=E+1; E1=E1+1;
                                end
                            end
                            
                        else
                            warndlg('No digital triggers to import on 2nd digital channel','Import .TRC Warning!');
                            for E=1:length(Trigs1)
                                TRC.event(E).type=PARAM.loadevents.dig_ch1_label;
                                TRC.event(E).latency=Trigs1(E);
                                TRC.event(E).urevent=E;
                                TRC.urevent(E).type=PARAM.loadevents.dig_ch1_label;
                                TRC.urevent(E).latency=Trigs1(E);
                            end
                        end
                    end
                end
        end
end

fprintf('\nPreparing output...\n');
if PARAM.chan_adjust_status==1 && strcmpi(PARAM.loaddata.state,'yes')
    if length(mmtrace)>(fs*60)
        avgvar=mean(var(mmtrace(1:Num_Chan-3,1:fs*60)'));
    else
        avgvar=mean(var(mmtrace(1:Num_Chan-3,:)'));
    end
    
    ch_adj_t=['[' PARAM.chan_adjust ']'];
    ch_adj=eval(ch_adj_t);
    
    if ~strcmpi('PARAM.loaddata.type','raw')
        for ch=1:length(ch_adj)
            mmtrace(ch_adj(ch),:)=mmtrace(ch_adj(ch),:)*avgvar/var(mmtrace(ch_adj(ch),:));
        end
    end
end

chan_del=0;
if PARAM.chan_exclude_status==1 && strcmpi(PARAM.loaddata.state,'yes')
    if strcmpi(PARAM.movenonEEGchannels,'yes')
        error('cannot delete channels if markerchannels are moved. Might be fixed some day.');
    end
    ch_exl_t=['[' PARAM.chan_exclude ']'];
    ch_exl=eval(ch_exl_t);
    e=1;c=1;shift=0;
    mmtrace(ch_exl,:)=[];
    chan_del=length(ch_exl);
end

if strcmpi(PARAM.loaddata.state,'yes')
    if strcmpi(PARAM.loaddata.type,'raw')
        
        if Bytes==2
            TRC.data=mmtrace-32768;
            disp('data is read as 16 bit raw data')
        else
            TRC.data=mmtrace-2097152;
            disp('data is read as 22 bit raw data')
        end
    else
        TRC.data=-mmtrace*1e6; % scale to uV and change polarity for EEGLAB
    end
end

sp=findstr(surname,'  ');
TRC.filename=char(PARAM.filename);
if sp >=1
    TRC.setname=[surname(1:(sp-1)) ', ' name(1) '. ' year month day ' .TRC File'];
else
    TRC.setname=[surname ', ' name(1) '. ' year month day ' .TRC File'];
end
TRC.starttime=[day '-' num2str(monthnr) '-' year ' ' num2str(hour,'%02.2d') ':' num2str(minutes,'%02.2d') ':' num2str(seconds,'%02.2d')];
TRC.time=m/fs;
timeins=3600*hour+60*minutes+seconds+fix(m/fs);
hour=fix(timeins/3600);
timeins=mod(timeins,3600);
minutes=fix(timeins/60);
seconds=mod(timeins,60);
TRC.endtime=[day '-' num2str(monthnr) '-' year ' ' num2str(hour,'%02.2d') ':' num2str(minutes,'%02.2d') ':' num2str(seconds,'%02.2d')];
TRC.pnts=m;
TRC.xmax=(TRC.pnts-1)/fs;
TRC.filename=trcfile;
TRC.filepath='';

TRC.nbchan=Num_Chan-chan_del;
TRC.trials=1;
TRC.srate=fs;
TRC.xmin=0;
TRC.ref='common';
TRC.bytes=Bytes;
TRC.channel=channel;
TRC.channelx=channelx;
for c=1:TRC.nbchan
    TRC.chanlocs(c).labels=electrode(c).label;
end

if strcmpi(PARAM.verbose,'yes')
    TRC.title=trctitle;
    TRC.laboratory=laboratory;
    TRC.birth_month=birth_month;
    TRC.birth_day=birth_day;
    TRC.birth_year=birth_year;
    TRC.acquisition_unit=acquisition_unit;
    TRC.acquisition_unit_description=acquisition_unit_description;
    TRC.filetype=filetype;
    TRC.filetype_code=filetype_code;
    TRC.filetype_description=filetype_description;
    TRC.Code_Area=Code_Area;
    TRC.Electrode_Area=Electrode_Area;
    TRC.electrode=electrode;
    TRC.originalelectrode=originalelectrode;
    
end






return;