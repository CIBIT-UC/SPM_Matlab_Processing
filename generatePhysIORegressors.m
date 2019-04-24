function [Regressors] = generatePhysIORegressors(file_path,file_prefix,data_configs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

physio = tapas_physio_new();

physio.log_files.vendor = 'Siemens_Tics';
physio.log_files.cardiac = {fullfile(file_path,['Physio_' file_prefix '_PULS.log'])};
% physio.log_files.respiration = {fullfile(file_path,['Physio_' file_prefix '_RESP.log'])};
physio.log_files.scan_timing = {fullfile(file_path,['Physio_' file_prefix '_Info.log'])};
physio.log_files.relative_start_acquisition = 0;
physio.log_files.align_scan = 'last';
physio.scan_timing.sqpar.Nslices = data_configs.nslices;
physio.scan_timing.sqpar.TR = data_configs.tr;
physio.scan_timing.sqpar.Ndummies = data_configs.ndummies;
physio.scan_timing.sqpar.Nscans = data_configs.nscans;
physio.scan_timing.sqpar.onset_slice = data_configs.onsetslice;
physio.scan_timing.sync.method = 'scan_timing_log';
physio.preproc.cardiac.modality = 'PPU';
physio.preproc.cardiac.initial_cpulse_select.method = 'auto_matched';
physio.preproc.cardiac.initial_cpulse_select.file = 'initial_cpulse_kRpeakfile.mat';
physio.preproc.cardiac.initial_cpulse_select.min = 0.4;
physio.preproc.cardiac.posthoc_cpulse_select.method = 'off';
physio.preproc.cardiac.posthoc_cpulse_select.percentile = 80;
physio.preproc.cardiac.posthoc_cpulse_select.upper_thresh = 60;
physio.preproc.cardiac.posthoc_cpulse_select.lower_thresh = 60;
physio.model.orthogonalise = 'none';
physio.model.censor_unreliable_recording_intervals = false;
% physio.model.output_multiple_regressors = 'multiple_regressors.txt';
% physio.model.output_physio = 'physio.mat';
physio.model.retroicor.include = true;
physio.model.retroicor.order.c = 3;
physio.model.retroicor.order.r = 4;
physio.model.retroicor.order.cr = 1;

% matlabbatch{1}.spm.tools.physio.model.rvt.yes.delays = 0;
% matlabbatch{1}.spm.tools.physio.model.hrv.yes.delays = 0;

physio.model.rvt.include = false;
physio.model.rvt.delays = 0;
physio.model.hrv.include = false;
physio.model.hrv.delays = 0;
physio.model.noise_rois.include = false;
physio.model.noise_rois.thresholds = 0.9;
physio.model.noise_rois.n_voxel_crop = 0;
physio.model.noise_rois.n_components = 1;
physio.model.movement.include = false;
physio.model.movement.order = 6;
physio.model.movement.censoring_threshold = 0.5;
physio.model.movement.censoring_method = 'FD';
physio.model.other.include = false;
physio.verbose.level = 2;
physio.verbose.process_log = cell(0, 1);
physio.verbose.fig_handles = zeros(0, 1);
physio.verbose.use_tabs = false;
physio.ons_secs.c_scaling = 1;
physio.ons_secs.r_scaling = 1;

[~, Regressors, ~] = tapas_physio_main_create_regressors(physio);

end
