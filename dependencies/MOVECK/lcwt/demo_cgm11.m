%% DEMO FOR THE CGM 1.1 model
%  Read the help of the function 'compute_cgm11' to get details of the
%  structures and fields used below

%% Trial filenames
% Set the path of the files to use for the calibration and the task to
% analyze respectively
trial_filenames.calib_mocap = 'C:\Users\SEAD\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\HUG2020\Programme\Test_Moveck_ComputationCGM\DATA_TEST\03080_04976_20190826-SBNNN-VDEF-02.C3D';
trial_filenames.task_mocap = 'C:\Users\SEAD\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\HUG2020\Programme\Test_Moveck_ComputationCGM\DATA_TEST\03080_04976_20190826-GBNNN-VDEF-08.C3D';
%% Sujects details
% Adapt the values below based on the anthropometry of your subject
subject_details.body_region = 'FullBody';
subject_details.marker_diameter = 16;
subject_details.subject_mass = 47.5;
subject_details.head_offset_enabled = true;
subject_details.left_shoulder_offset = 50;
subject_details.right_shoulder_offset = 50;
subject_details.left_elbow_width = 80;
subject_details.right_elbow_width = 80;
subject_details.left_wrist_width = 40;
subject_details.right_wrist_width = 40;
subject_details.left_hand_thickness = 30;
subject_details.right_hand_thickness = 30;
subject_details.left_leg_length = 800;
subject_details.right_leg_length = 800;
subject_details.left_knee_width = 105;
subject_details.right_knee_width = 105;
subject_details.left_ankle_width = 55;
subject_details.right_ankle_width = 55;
subject_details.left_foot_flat_enabled = true;
subject_details.right_foot_flat_enabled = true;
% Adapt the values below to assign the forceplate with the feet
subject_details.forceplate_index_left_foot_map = 2; 
subject_details.forceplate_index_right_foot_map = 1;
% CGM1.1 computation 
tstart = tic;
[output, store] = compute_cgm11(trial_filenames, subject_details);
disp(['Total time: ', num2str(toc(tstart))]);
%% Data comparisons
if(0)
compare_subject_center_of_mass(store, output)
compare_subject_joint_angles(store, output)
compare_subject_joint_force(store, output)
compare_subject_joint_moment(store, output)
compare_subject_grw(store, output); % must be executed after force and moment
end
%% Data statistics
plot_kinematics_cycles_lowerlimb(store, output);
%%
