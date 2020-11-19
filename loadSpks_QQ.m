function [dataRaw, dataClean, dataOutlier] = loadSpks_QQ(i_name)
%loadSpks_QQ: load spikes from R. Q. Quiroga, Z. Nadasdy, and Y. Ben-Shaul
% Unsupervised spike detection and sorting with wavelets and superparamagnetic clustering 
% Neural Comput., vol. 16, no. 8, pp. 1661–1687, 2004.
%
% dataRaw/dataClean
%   spikes: m_variable * n_observation
%   spike_class: 1 * n_obseravation

datafile_names = {...
    'C_Easy1_noise005', 'C_Easy1_noise010', 'C_Easy1_noise015', 'C_Easy1_noise020',  'C_Easy1_noise025', 'C_Easy1_noise030', 'C_Easy1_noise035', 'C_Easy1_noise040', ...
    'C_Easy2_noise005', 'C_Easy2_noise010', 'C_Easy2_noise015', 'C_Easy2_noise020',...
    'C_Difficult1_noise005', 'C_Difficult1_noise010', 'C_Difficult1_noise015', 'C_Difficult1_noise020',...
    'C_Difficult2_noise005', 'C_Difficult2_noise010', 'C_Difficult2_noise015', 'C_Difficult2_noise020'
};
load( ['.\',datafile_names{i_name},'.mat'] );
Data.name = datafile_names{i_name};

Data.frequency = 1000/samplingInterval; % samplingInterval -- 1000ms/24000Hz
Data.data = data;
Data.OVERLAP_DATA = OVERLAP_DATA;
Data.spike_class = spike_class;
Data.spike_times = spike_times{1};
Data.chan = chan;
Data.startData = startData;

% filtering
Data.data = spikefilter(data);
%whether has it or not will not affect the results significantly.

dim_spk = 64;
spk_begin = Data.startData + Data.spike_times;
spk_end = spk_begin+dim_spk-1;
spks = zeros( dim_spk, numel(spk_begin) ); % 64_variable * n_observation
for i = 1:numel(spk_begin)
	spks(:,i) = Data.data( spk_begin(i) : spk_end(i) );
end
dataRaw.spikes = spks; %dataRaw
dataRaw.spike_class = Data.spike_class{1};
nonOverlap_indexes = find(Data.spike_class{2}==0); %dataClean
olpIdx = setdiff(1:length(spks),nonOverlap_indexes);
dataClean.spikes = spks(:,nonOverlap_indexes);
dataClean.spike_class = dataRaw.spike_class(nonOverlap_indexes);
dataOutlier.spike_class = Data.spike_class{1}; %dataOutlier
dataOutlier.spike_class(Data.spike_class{2}==1) = 0;
dataOutlier.spikes = spks; 
fprintf( '===> Data %-2d: %-23s \n', i_name, Data.name );
end

function xf_detect = spikefilter(x, par)
% x: 1*1440000 raw data
% par
%   par.sr=24000; ampling rate (in Hz).
%   par.detect_fmin=300; high pass filter for sorting
%   par.detect_fmax=3000; low pass filter for sorting
%   par.detect_order=4;  filter order for detection (2 for extract spikes)
% ----
% xf_detect = 1*144000
% 
% Chaure, Fernando J., Hernan G. Rey, and Rodrigo Quian Quiroga.
% "A novel and fully automatic spike-sorting implementation with variable number of features."
% Journal of neurophysiology 120.4 (2018): 1859-1871.
if nargin<2
    par.sr = 24000;
    par.detect_fmin = 300;
    par.detect_fmax = 3000;
    par.detect_order = 2;
end    
    
sr = par.sr;
fmin_detect = par.detect_fmin;
fmax_detect = par.detect_fmax;
detect_order = par.detect_order;


% HIGH-PASS FILTER OF THE DATA
if exist('ellip','file') %Checks for the signal processing toolbox
    [b,a] = ellip(detect_order,0.1,40,[fmin_detect fmax_detect]*2/sr);
    if exist('FiltFiltM','file')
        xf_detect = FiltFiltM(b, a, x);
    else
        xf_detect = filtfilt(b, a, x);
    end
else
    xf_detect = fix_filter(x); %Does a bandpass filtering between [300 3000] without the toolbox.
end
end