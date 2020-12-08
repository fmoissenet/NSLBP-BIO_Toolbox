% Compute the joint kinematics based on the CGM 1.1 definition.
%
% The transformation steps are listed below:
%  - Load the data
%    - Markers trajectory for the calibration
%    - Markers trajectory for the analyzed task
%  - Prepare the data
%    - Data processing for the markers (see the function prepare_mocap)
%    - Data processing for the forceplate (see the function prepare_forceplate)
%  - Model calibration
%  - Model pose estimation
%  - Model joint kinematics estimation
%  - Matlab export (joint kinematics)
%  - Model relative BSiP
%  - Model global BSiP
%  - Model body kinematics
%  - Ground segment construction (from force platforms)
%  - Model body dynamics
%  - Model joint kinetics
%
% INPUTS
%  - trial_filenames: structure with the filename of the trials to load
%   - trial_filenames.calib_mocap: path to a C3D file used to calibrate the model
%   - trial_filenames.task_mocap: path to a C3D file containing the markers' coordinates of the task to analyze
%  - subject_details: structure with subject anthropometry
%   - body_region: string 'FullBody', 'LowerBody', 'UpperBody' to determine the part of the body to analyze
%   - marker_diameter: diameter of the marker used in the data acquisition (unit: mm)
%   - subject_mass: mass of the subject (in kg)
%   - head_offset_enabled: boolean to compute the offset angle of the head
%   - left_shoulder_offset: distance to compute the left shoulder joint centre (unit: mm)
%   - right_shoulder_offset: distance to compute the right shoulder joint centre (unit: mm)
%   - left_elbow_width: width of the left elbow (unit: mm)
%   - right_elbow_width: width of the right elbow (unit: mm)
%   - left_wrist_width: width of the left wrist (unit: mm)
%   - right_wrist_width: width of the right wrist (unit: mm)
%   - left_hand_thickness: thickness of the left hand (unit: mm)
%   - right_hand_thickness: thickness of the right hand (unit: mm)
%   - left_leg_length: length of the left leg (unit: mm)
%   - right_leg_length: length of the right leg (unit: mm)
%   - left_knee_width: width of the left knee (unit: mm)
%   - right_knee_width: width of the right knee (unit: mm)
%   - left_ankle_width: width of the left ankle (unit: mm)
%   - right_ankle_width: (unit: mm)
%   - left_foot_flat_enabled: boolean to set the left foot as flat
%   - right_foot_flat_enabled: boolean to set the right foot as flat
%
% OUTPUTS
%  - output: structure with result of the model
%    - output.trial_filenames: copy of the input 'trial_filenames'
%    - output.subject_details: copy of the input 'subject_details'
%    - output.task_analysis: structure with the joint kinematics associated with the analyzed task
%      - output.task_analysis.(Left|Right)HeadProgressionAngles: clinical angles for the head segment
%      - output.task_analysis.(Left|Right)NeckAngles: clinical angles for the shoulder joint
%      - output.task_analysis.(Left|Right)ShoulderAngles: clinical angles for the shoulder joint
%      - output.task_analysis.(Left|Right)ElbowAngles: clinical angles for the elbow joint
%      - output.task_analysis.(Left|Right)WristAngles: clinical angles for the wrist joint
%      - output.task_analysis.(Left|Right)TorsoProgressionAngles: clinical angles for the torso segment
%      - output.task_analysis.(Left|Right)SpineAngles: clinical angles for the shoulder joint
%      - output.task_analysis.(Left|Right)PelvisProgressionAngles: clinical angles for the pelvis segment
%      - output.task_analysis.(Left|Right)HipAngles: clinical angles for the hip joint
%      - output.task_analysis.(Left|Right)KneeAngles: clinical angles for the knee joint
%      - output.task_analysis.(Left|Right)AnkleAngles: clinical angles for the ankle joint
%      - output.task_analysis.(Left|Right)FootProgressionAngles: clinical angles for the (left|right) foot segment
%      - output.task_analysis.(Left|Right)HipForce: clinical force for the hip joint
%      - output.task_analysis.(Left|Right)KneeForce: clinical force for the knee joint
%      - output.task_analysis.(Left|Right)AnkleForce: clinical force for the ankle joint
%      - output.task_analysis.(Left|Right)HipMoment: clinical moment for the hip joint
%      - output.task_analysis.(Left|Right)KneeMoment: clinical moment for the knee joint
%      - output.task_analysis.(Left|Right)AnkleMoment: clinical moment for the ankle joint
%  - store: internal structure containing all the results (can be used for debugging or validation) 

% HISTORY
% Version 2.1 (2020-06-24)
%  - Progression angles added
%  - Neck angles added (for UpperBody and FullBody configuration)
%  - Spine angles added (for FullBody configuration only)
% Version 2.0 (2020-06-01)
%  - Computation of the joint kinetics and all required steps added
%  - Name for the groups/sets/attributes modified to follow the updated naming convention
% Version: 1.0 (2020-04-03)
%  - Initial version 
%
% Author: Arnaud Barre (Moveck Solution inc.) <arnaud.barre@moveck.com>
function [output, store] = compute_cgm11(trial_filenames, subject_details)
% Data loading
tstart=tic;
store = moveck.data_store();
disp([' * storage construction (might include plugins loading): ', num2str(toc(tstart))]);
tstart=tic;
trial_calib = moveck.import_trial(store, trial_filenames.calib_mocap);
disp([' * loading calibration file: ', num2str(toc(tstart))]);
tstart=tic;
trial_task = moveck.import_trial(store, trial_filenames.task_mocap);
disp([' * loading task file: ', num2str(toc(tstart))]);
% Data preparation
disp(' * mocap preparation for the calibration file:');
tstart=tic;
prepare_mocap(trial_calib);
disp(['    > TOTAL: ', num2str(toc(tstart))]);
disp(' * mocap preparation for the task file:');
tstart=tic;
prepare_mocap(trial_task);
disp(['    > TOTAL: ', num2str(toc(tstart))]);
disp(' * forceplate preparation for the task file:');
tstart=tic;
prepare_forceplate(trial_task, trial_task);
disp(['    > TOTAL: ', num2str(toc(tstart))]);
destination_group = 'Models/CGM11';
landmarks_map = {...
  'LFHD', 'LHF', ...
  'LBHD', 'LHB', ...
  'RFHD', 'RHF', ...
  'RBHD', 'RHB', ...
  'C7'  , 'C7', ...  % same
  'T10' , 'T10', ... % same
  'CLAV', 'SS', ...
  'STRN', 'XP', ...
  'LSHO', 'LAC', ...
  'RSHO', 'RAC', ...
  'LELB', 'LLHE', ...
  'LWRB', 'LUS', ...
  'LWRA', 'LRS', ...
  'RELB', 'RLHE', ...
  'RWRB', 'RUS', ...
  'RWRA', 'RRS', ...
  'LFIN', 'LMH2', ...
  'RFIN', 'RMH2', ...
  'LASI', 'LASIS', ...
  'RASI', 'RASIS', ...
  'LPSI', 'LPSIS', ...
  'RPSI', 'RPSIS', ...
  'SACR', 'SC', ...
  'LTHI', 'LITB', ...
  'RTHI', 'RITB', ...
  'LKNE', 'LLFE', ...
  'RKNE', 'RLFE', ...
  'LTIB', 'LLS', ...
  'RTIB', 'RLS', ...
  'LANK', 'LLTM', ...
  'RANK', 'RLTM', ...
  'LTOE', 'LMTH2', ...
  'RTOE', 'RMTH2', ...
  'LHEE', 'LHEE', ... % same
  'RHEE', 'RHEE'};    % same
% Model calibration
static_calibration_settings = struct();
static_calibration_settings.BodyRegion = subject_details.body_region;
static_calibration_settings.LandmarksMap = landmarks_map;
static_calibration_settings.SourceGroup = 'Processings/Mocap';
static_calibration_settings.DestinationGroup = destination_group;
static_calibration_settings.MarkerDiameter = subject_details.marker_diameter;
static_calibration_settings.HeadOffsetEnabled = subject_details.head_offset_enabled;
static_calibration_settings.LeftShoulderOffset = subject_details.left_shoulder_offset;
static_calibration_settings.RightShoulderOffset = subject_details.right_shoulder_offset;
static_calibration_settings.LeftElbowWidth = subject_details.left_elbow_width;
static_calibration_settings.RightElbowWidth = subject_details.right_elbow_width;
static_calibration_settings.LeftWristWidth = subject_details.left_wrist_width;
static_calibration_settings.RightWristWidth = subject_details.right_wrist_width;
static_calibration_settings.LeftHandThickness = subject_details.left_hand_thickness;
static_calibration_settings.RightHandThickness = subject_details.right_hand_thickness;
static_calibration_settings.LeftLegLength = subject_details.left_leg_length;
static_calibration_settings.RightLegLength = subject_details.right_leg_length;
static_calibration_settings.LeftKneeWidth = subject_details.left_knee_width;
static_calibration_settings.LeftAnkleWidth = subject_details.left_ankle_width;
static_calibration_settings.LeftFootFlatEnabled = subject_details.left_foot_flat_enabled;
static_calibration_settings.RightKneeWidth = subject_details.right_knee_width;
static_calibration_settings.RightAnkleWidth = subject_details.right_ankle_width;
static_calibration_settings.RightFootFlatEnabled = subject_details.right_foot_flat_enabled;
if (isfield(subject_details, 'inter_asis_distance'))
    static_calibration_settings.InterAsisDistance = subject_details.inter_asis_distance;
end
static_calibration_settings.callable_unit = 'cgm1x.calibration';
tstart=tic;
moveck.transform_data(trial_calib, static_calibration_settings);
disp([' * CGM 1.1 calibration: ', num2str(toc(tstart))]);
% Copy the model calibration details to the task trial
copy_group_settings = struct();
copy_group_settings.SourceGroup = [trial_calib.name(), '/', static_calibration_settings.DestinationGroup];
copy_group_settings.DestinationGroup = static_calibration_settings.DestinationGroup;
copy_group_settings.ScanRecursive = true;
copy_group_settings.callable_unit = 'data-modifier.group-copy';
tstart=tic;
moveck.transform_data(trial_task, copy_group_settings);
disp([' * CGM 1.1 group copy: ', num2str(toc(tstart))]);
% Model pose estimation
pose_estimation_settings = struct();
pose_estimation_settings.LandmarksMap = landmarks_map;
pose_estimation_settings.SourceGroup = 'Processings/Mocap';
pose_estimation_settings.DestinationGroup = destination_group;
pose_estimation_settings.callable_unit = 'cgm1x.reconstruction';
tstart=tic;
moveck.transform_data(trial_task, pose_estimation_settings);
disp([' * CGM 1.1 reconstruction: ', num2str(toc(tstart))]);
% Joint kinematics estimation
joint_kinematics_settings = struct();
joint_kinematics_settings.SourceGroup = destination_group;
joint_kinematics_settings.DestinationGroup = destination_group;
joint_kinematics_settings.callable_unit = 'cgm1x.joint-angles';
tstart=tic;
moveck.transform_data(trial_task, joint_kinematics_settings);
disp([' * CGM 1.1 joint kinematics: ', num2str(toc(tstart))]);
% Progression axis detection
progression_axis_settings = struct();
progression_axis_settings.SourceGroup = destination_group;
progression_axis_settings.DestinationGroup = destination_group;
progression_axis_settings.SegmentReferenceHints = {'Pelvis', 'Torso'};
progression_axis_settings.callable_unit = 'classical-mechanics.straight-progression-axis';
tstart=tic;
moveck.transform_data(trial_task, progression_axis_settings);
disp([' * CGM 1.1 progression axis: ', num2str(toc(tstart))]);
% Progression kinematics estimation
progression_kinematics_settings = struct();
progression_kinematics_settings.SourceGroup = destination_group;
progression_kinematics_settings.DestinationGroup = destination_group;
progression_kinematics_settings.callable_unit = 'cgm1x.progression-angles';
tstart=tic;
moveck.transform_data(trial_task, progression_kinematics_settings);
disp([' * CGM 1.1 progression kinematics: ', num2str(toc(tstart))]);
% Export kinematics results in Matlab
tstart=tic;
output = struct();
output.trial_filenames = trial_filenames;
output.subject_details = subject_details;
outputs_model = trial_task.retrieve_group(destination_group);
joints_name = outputs_model.retrieve_group('Joints').list_group_children_name();
for i = 1:length(joints_name)
    joint_kinematics_name = ['Joints/', joints_name{i}, '/Angles'];
    output.task_analysis.([joints_name{i}, 'Angles']) = squeeze(outputs_model.retrieve_set(joint_kinematics_name).read())';
end
disp([' * Matlab kinematics export: ', num2str(toc(tstart))]);
% Model relative BSiP
relative_body_inertial_settings = struct();
relative_body_inertial_settings.SourceGroup = destination_group;
relative_body_inertial_settings.DestinationGroup = destination_group;
relative_body_inertial_settings.SubjectMass = subject_details.subject_mass;
relative_body_inertial_settings.callable_unit = 'cgm1x.relative-inertial-parameters';
tstart=tic;
moveck.transform_data(trial_task, relative_body_inertial_settings);
disp([' * CGM 1.1 relative BSiP: ', num2str(toc(tstart))]);
% Model global BSiP
body_inertial_settings = struct();
body_inertial_settings.SourceGroup = destination_group;
body_inertial_settings.DestinationGroup = destination_group;
body_inertial_settings.callable_unit = 'classical-mechanics.inertial-parameters';
tstart=tic;
moveck.transform_data(trial_task, body_inertial_settings);
% FIXME: Wrong CoM definition for the torso. We have to update the
% Segments/Pelvis/p_COM coordinates
torso = trial_task.retrieve_group([destination_group, '/Segments/Torso']);
torso_prox = torso.retrieve_set('p_Proximal').read(); % C7 vertebra
torso_dist = trial_task.retrieve_set([destination_group, '/Segments/Pelvis/p_Proximal']).read(); % L5 vertebra
torso_rel_com = torso.retrieve_attribute('p_COM^SCS').read() / torso.retrieve_attribute('Length').read();
torso_com = torso_prox + (torso_dist - torso_prox) * torso_rel_com(3);
torso.retrieve_set('p_COM').write(torso_com);
disp([' * CGM 1.1 global BSiP: ', num2str(toc(tstart))]);
% Model whole body CoM
whole_centre_mass_settings = struct();
whole_centre_mass_settings.SourceGroup = destination_group;
whole_centre_mass_settings.DestinationGroup = destination_group;
whole_centre_mass_settings.callable_unit = 'classical-mechanics.whole-centre-mass';
tstart=tic;
moveck.transform_data(trial_task, whole_centre_mass_settings);
disp([' * CGM 1.1 Whole body CoM: ', num2str(toc(tstart))]);
% Model body kinematics
body_kinematics_settings = struct();
body_kinematics_settings.SourceGroup = destination_group;
body_kinematics_settings.DestinationGroup = destination_group;
body_kinematics_settings.callable_unit = 'classical-mechanics.body-kinematics';
tstart=tic;
moveck.transform_data(trial_task, body_kinematics_settings);
disp([' * CGM 1.1 body kinematics: ', num2str(toc(tstart))]);
% % Ground segments
% tstart=tic;
% kinetics_chains_processed = {};
% if (isfield(subject_details, 'forceplate_index_left_foot_map'))
%     left_wrench = trial_task.retrieve_set(['Processings/FP',num2str(subject_details.forceplate_index_left_foot_map),'/W_POA']).read();
%     left_ground_group = trial_task.create_group([destination_group,'/Segments/LeftGround']);
%     left_ground_group.create_set('F_Proximal', -left_wrench(1:3,:,:));
%     left_ground_group.create_set('M_Proximal', -left_wrench(4:6,:,:));
%     left_ground_group.create_set('p_Proximal', left_wrench(7:9,:,:));
%     kinetics_chains_processed{end+1} = 'LeftLowerLimb';
% end
% if (isfield(subject_details, 'forceplate_index_right_foot_map'))
%     right_wrench = trial_task.retrieve_set(['Processings/FP',num2str(subject_details.forceplate_index_right_foot_map),'/W_POA']).read();
%     right_ground_group = trial_task.create_group([destination_group,'/Segments/RightGround']);
%     right_ground_group.create_set('F_Proximal', -right_wrench(1:3,:,:));
%     right_ground_group.create_set('M_Proximal', -right_wrench(4:6,:,:));
%     right_ground_group.create_set('p_Proximal', right_wrench(7:9,:,:));
%     kinetics_chains_processed{end+1} = 'RightLowerLimb';
% end
% disp([' * Ground reaction wrench(es) assignment: ', num2str(toc(tstart))]);
% % Body dynamics
% body_dynamics_settings = struct();
% body_dynamics_settings.SourceGroup = destination_group;
% body_dynamics_settings.DestinationGroup = destination_group;
% body_dynamics_settings.GravityDirection = [0,0,-1];
% body_dynamics_settings.Chains = kinetics_chains_processed;
% body_dynamics_settings.callable_unit = 'classical-mechanics.body-dynamics';
% tstart=tic;
% moveck.transform_data(trial_task, body_dynamics_settings);
% disp([' * CGM 1.1 body dynamics: ', num2str(toc(tstart))]);
% % Joint kinetics
% tstart=tic;
% markers_num_samples = trial_task.retrieve_group('Format/Config/Points').retrieve_attribute('NumSamples').read();
% wrench_names = {'LeftFoot', 'LeftShank', 'LeftThigh', 'RightFoot', 'RightShank', 'RightThigh'};
% segment_names = {'LeftShank', 'LeftThigh', 'Pelvis', 'RightShank', 'RightThigh', 'Pelvis'};
% joint_names = {'LeftAnkle', 'LeftKnee', 'LeftHip', 'RightAnkle', 'RightKnee', 'RightHip'};
% interpretation_force = {{3,2,-1}, {1,2,3}, {-1,2,3}, {3,2,-1}, {1,2,3}, {-1,2,3}};
% interpretation_moment = {{2,1,3}, {2,-1,3}, {2,1,-3}, {2,1,-3}, {-2,-1,-3}, {-2,1,-3}};
% for i = 1:length(wrench_names)
%     wrench_name = wrench_names{i};
%     segment_name = segment_names{i};
%     joint_name = joint_names{i};
%     f_name = ['Segments/', wrench_name, '/F_Proximal'];
%     m_name = ['Segments/', wrench_name, '/M_Proximal'];
%     scs_name = ['Segments/', segment_name, '/T_SCS'];
%     if (~outputs_model.exists_set(f_name) || ~outputs_model.exists_set(m_name) || ~outputs_model.exists_set(scs_name))
%         fprintf('Missing set %s, %s, or %s\n', f_name, m_name, scs_name);
%         continue;
%     end
%     scs = outputs_model.retrieve_set(scs_name).read();
%     f = outputs_model.retrieve_set(f_name).read();
%     m = outputs_model.retrieve_set(m_name).read();
%     if (size(scs, 3) ~= markers_num_samples || size(scs, 3) ~= size(f, 3) || size(scs, 3) ~= size(m, 3))
%         error('lcwt:cgm11', ['Incompatible number of samples to compute joint kinetics for the joint ''', joint_name, '''']);
%     end
%     forces = nan * zeros(markers_num_samples, 3);
%     moments = nan * zeros(markers_num_samples, 3);
%     interp_force = interpretation_force{i};
%     interp_moment = interpretation_moment{i};
%     for j = 1:markers_num_samples
%         R = transpose(scs(1:3, 1:3, j));
%         f_ = R * f(:, 1, j);
%         m_ = R * m(:, 1, j);
%         for k = 1:3
%             idx = abs(interp_force{k});
%             s = sign(interp_force{k});
%             forces(j, idx) = s * f_(k);
%             idx = abs(interp_moment{k});
%             s = sign(interp_moment{k});
%             moments(j, idx) = s * m_(k);
%         end
%     end
%     output.task_analysis.([joint_name, 'Force']) = forces;
%     output.task_analysis.([joint_name, 'Moment']) = moments;
% end
% disp([' * CGM 1.1 joint kinetics: ', num2str(toc(tstart))]);
% Events extraction (available in a C3D file)
% c3d_events_detect_settings = struct();
% c3d_events_detect_settings.DestinationGroup = 'Analyses/CGM11/Events';
% c3d_events_detect_settings.callable_unit = 'org.c3d.events-detect';
% tstart=tic;
% moveck.transform_data(trial_task, c3d_events_detect_settings);
% disp([' * Events extraction: ', num2str(toc(tstart))]);
% % Cycles indexation and parameters preparation
% index_from_events_settings = struct();
% index_from_events_settings.SourceGroup = 'Analyses/CGM11';
% index_from_events_settings.DestinationGroup = 'Analyses/CGM11';
% index_from_events_settings.callable_unit = 'cycle.index_from_events';
% parameters_preparation_settings = struct();
% parameters_preparation_settings.SourceGroup = 'Models/CGM11';
% parameters_preparation_settings.DestinationGroup = 'Analyses/CGM11';
% parameters_preparation_settings.callable_unit = 'cycle.prepare_parameters';
% % - Left
% index_from_events_settings.Events = {'LHS'};%{'LeftFootStrike'};
% index_from_events_settings.Context = 'General';%'Left';
% tstart=tic;
% moveck.transform_data(trial_task, index_from_events_settings);
% disp([' * Left cycle indexation: ', num2str(toc(tstart))]);
% parameters_preparation_settings.Context = 'Left';
% parameters_preparation_settings.ParametersMap = {...
%     'Joints/LeftAnkle/Angles', 'LeftAnkleAngles', ...
%     'Joints/LeftKnee/Angles', 'LeftKneeAngles', ...
%     'Joints/LeftHip/Angles', 'LeftHipAngles', ...
% };
% tstart=tic;
% moveck.transform_data(trial_task, parameters_preparation_settings);
% disp([' * Left parameters preparation: ', num2str(toc(tstart))]);
% % - Right
% index_from_events_settings.Events = {'RHS'};%{'RightFootStrike'};
% index_from_events_settings.Context = 'General';%'Right';
% tstart=tic;
% moveck.transform_data(trial_task, index_from_events_settings);
% disp([' * Right cycle indexation: ', num2str(toc(tstart))]);
% parameters_preparation_settings.Context = 'Right';
% parameters_preparation_settings.ParametersMap = {...
%     'Joints/RightAnkle/Angles', 'RightAnkleAngles', ...
%     'Joints/RightKnee/Angles', 'RightKneeAngles', ...
%     'Joints/RightHip/Angles', 'RightHipAngles', ...
% };
% tstart=tic;
% moveck.transform_data(trial_task, parameters_preparation_settings);
% disp([' * Right parameters preparation: ', num2str(toc(tstart))]);
% % Time linear length normalization
% time_linear_length_normalization_settings = struct();
% time_linear_length_normalization_settings.SourceGroup = 'Analyses/CGM11';
% time_linear_length_normalization_settings.DestinationGroup = 'Analyses/CGM11';
% time_linear_length_normalization_settings.callable_unit = 'cycle.linear_length_normalization';
% tstart=tic;
% moveck.transform_data(trial_task, time_linear_length_normalization_settings);
% disp([' * Linear length normalization: ', num2str(toc(tstart))]);