% Prepare the content of a trial to compute a ground reaction wrench
%
% The transformation steps are listed below:
%  - Split the original set to channels
%  - Copy the channels to the mocap trial
%  - Downsample the content of each copied channel
%  - Filter the channels with a zero-lag butterworth low-pass filter
%  - Create a force plateform
%  - Compute the associated ground reaction wrench 
%
% INPUTS
%  - trial_forceplate: the original trial with the raw forceplate signals
%  - trial_mocap: (optional) the destination trial with mocap data

% Version: 1.5 (2020-06-11)
%  - Callable 'set-filter-zero-lag-butterworth' replaced by 'signal.set-butterworth'
% Version: 1.4 (2020-05-26)
%  - The forceplates detection and GRWs computation are enabled
%  - The Butterworth low-pass filter is done before the downsampling to
%    improve the quality (smoothness) of the signal
%  - Specifications' names updated to use the CamelCase convention 
% Version: 1.3 (2020-04-01)
%  - Custom functions (e.g. 'moveck.split-set') was replaced by the standard
%    way using the identifier of the callable unit to use and the function
%    'moveck.transform_data'
%  - The forceplates detection and GRWs computation are disabled temporary
% Version: 1.2 (2020-03-01)
%  - Initial version 
%
% Author: Arnaud Barr� (Moveck Solution inc.)
function prepare_forceplate(trial_forceplate, trial_mocap)
    if (nargin == 1)
        trial_mocap = trial_forceplate;
    end
    %% I. Prepare data
    % I.1. Create the 'adc' instrument (split the tensor 'analogs')
    analogs_split_settings = struct();
    analogs_split_settings.SourceSet = 'Format/Data/Analogs';
    analogs_split_settings.DestinationGroup = 'Instruments/ADC';
    analogs_split_settings.callable_unit = 'data-modifier.set-split';
    tstart=tic;
    moveck.transform_data(trial_forceplate, analogs_split_settings);
    disp(['    > data split to create ADC: ', num2str(toc(tstart))]);
    % I.2. Copy the content of the 'adc' instrument to be processed
    copy_group_settings = struct();
    copy_group_settings.SourceGroup = [trial_forceplate.name(), '/', analogs_split_settings.DestinationGroup];
    copy_group_settings.DestinationGroup = 'Processings/ADC';
    copy_group_settings.callable_unit = 'data-modifier.group-copy';
    tstart=tic;
    moveck.transform_data(trial_mocap, copy_group_settings);
    disp(['    > ADC data copied to be processed: ', num2str(toc(tstart))]);
    %% II. Signal processing
    % II.1. Zero-lag Butterworth low-pass filter
    zlbw_settings = struct();
    zlbw_settings.SourceGroup = copy_group_settings.DestinationGroup;
    zlbw_settings.CutoffFrequency = 20;
    zlbw_settings.Order = 4;
    zlbw_settings.Bandform = 'ZeroLagLowPass';
    zlbw_settings.callable_unit = 'signal.set-butterworth';
    tstart = tic;
    moveck.transform_data(trial_mocap, zlbw_settings);
    disp(['    > ADC zero lag low-pass buterworth filter: ', num2str(toc(tstart))]);
    % II.1. Downsample 
    downsample_settings = struct();
    downsample_settings.SourceGroup = copy_group_settings.DestinationGroup;
    downsample_settings.HintGroup = 'Processings/Mocap';
    downsample_settings.callable_unit = 'signal.set-downsample';
    tstart = tic;
    moveck.transform_data(trial_mocap, downsample_settings);
    disp(['    > ADC downsample: ', num2str(toc(tstart))]);
    %% III. Detect the C3D force plate(s) and use the processed signals as inputs
    fp_settings = struct();
    fp_settings.callable_unit = 'org.c3d.forceplate-detect';
    fp_settings.SourceGroup =  [trial_mocap.name(), '/', copy_group_settings.DestinationGroup];
    % - FP1
    fp_settings.DestinationGroup = [trial_mocap.name(), '/', 'Processings/FP1'];
    fp_settings.ForceplateIndex = 0; % WARNING: 0-based index!
    tstart = tic;
    moveck.transform_data(trial_forceplate, fp_settings);
    disp(['    > FP1 detection: ', num2str(toc(tstart))]);
    % - FP2
    fp_settings.DestinationGroup = [trial_mocap.name(), '/', 'Processings/FP2'];
    fp_settings.ForceplateIndex = 1; % WARNING: 0-based index!
    tstart = tic;
    moveck.transform_data(trial_forceplate, fp_settings);
    disp(['    > FP2 detection: ', num2str(toc(tstart))]);
    %% IV. Compute the ground reaction wrench at the point of application
    %      - Fz threshold set by default to 10N
    %      - Expressed in the global frame
    grw_settings = struct();
    grw_settings.callable_unit = 'classical-mechanics.ground-reaction-wrench';
    grw_settings.Location = 'point-of-application';
    % - FP1
    grw_settings.SourceGroup = 'Processings/FP1';
    grw_settings.DestinationGroup = 'Processings/FP1';
    tstart = tic;
    moveck.transform_data(trial_mocap, grw_settings);
    disp(['    > FP1 GRW computation: ', num2str(toc(tstart))]);
    % - FP2
    grw_settings.SourceGroup = 'Processings/FP2';
    grw_settings.DestinationGroup = 'Processings/FP2';
    tstart = tic;
    moveck.transform_data(trial_mocap, grw_settings);
    disp(['    > FP1 GRW computation: ', num2str(toc(tstart))]);
end