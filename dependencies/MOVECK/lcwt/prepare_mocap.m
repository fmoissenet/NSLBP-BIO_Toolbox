% Prepare the content of a trial to get filtered marker trajectories
%
% The transformation steps are listed below:
%  - Split the original set to get markers trajectories
%  - Copy the concent to be processed
%  - Fill gap the markers
%  - Filter the markers with a zero-lag butterworth low-pass filter
%
% INPUTS
%  - trial: a trial with mocap data

% Version: 1.5 (2020-06-11)
%  - Callable 'set-fill-gap' replaced by 'signal.set-fill-gap'
%  - Callable 'set-filter-zero-lag-butterworth' replaced by 'signal.set-butterworth'
% Version: 1.4 (2020-06-03)
%  - Specifications' names updated to use the CamelCase convention
%  - Callable 'group-copy' replaced by 'data-modifier.group-copy'
% Version: 1.3 (2020-04-01)
%  - Custom functions (e.g. 'moveck.split-set') was replaced by the standard
%    way using the identifier of the callable unit to use and the function
%    'moveck.transform_data'
% Version: 1.2 (2020-03-13)
%  - Initial version 
%
% Author: Arnaud Barr� (Moveck Solution inc.)
function prepare_mocap(trial)
    %% I. Prepare data
    % I.1. Create the 'mocap' instrument
    %      (Split the tensor 'points' and keep only the sets known as 'marker')
    markers_split_settings = struct();
    markers_split_settings.SourceSet = 'Format/Data/Points';
    markers_split_settings.DestinationGroup = 'Instruments/Mocap';
    markers_split_settings.Filter = {'Types', 'marker'};
    markers_split_settings.callable_unit = 'data-modifier.set-split';
    tstart=tic;
    moveck.transform_data(trial, markers_split_settings);
    disp(['    > data split to create markers: ', num2str(toc(tstart))]);
    % 2. Copy the content of the 'mocap' instrument to be processed
    copy_group_settings = struct();
    copy_group_settings.SourceGroup = markers_split_settings.DestinationGroup;
    copy_group_settings.DestinationGroup = 'Processings/Mocap';
    copy_group_settings.callable_unit = 'data-modifier.group-copy';
    tstart=tic;
    moveck.transform_data(trial, copy_group_settings);
    disp(['    > markers data copied to be processed: ', num2str(toc(tstart))]);
    %% II. Signal processing
    % II.1. Fill gap
    fill_gap_settings = struct();
    fill_gap_settings.SourceGroup = copy_group_settings.DestinationGroup;
    fill_gap_settings.HintMaxGapLength = 0.1; %Factor between 0 and 1 that select the max gap len based on the sample frequency
    fill_gap_settings.callable_unit = 'signal.set-fill-gap';
    tstart=tic;
    moveck.transform_data(trial, fill_gap_settings);
    disp(['    > markers fill gap: ', num2str(toc(tstart))]);
    % II.2. Zero-lag Butterworth low-pass Ffilter
    zlbw_settings = struct();
    zlbw_settings.SourceGroup = copy_group_settings.DestinationGroup;
    zlbw_settings.CutoffFrequency = 6;
    zlbw_settings.Order = 4;
    zlbw_settings.Bandform = 'ZeroLagLowPass';
    zlbw_settings.callable_unit = 'signal.set-butterworth';
    tstart=tic;
    moveck.transform_data(trial, zlbw_settings);
    disp(['    > markers zero lag low-pass buterworth filter: ', num2str(toc(tstart))]);
end