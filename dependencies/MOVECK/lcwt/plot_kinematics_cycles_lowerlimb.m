function plot_kinematics_cycles_lowerlimb(store, output)
    side_names = {'Left', 'Right'};
    lowerlimb_angle_names = {'HipAngles', 'KneeAngles', 'AnkleAngles'};
    angle_suffixes = {'F/E', 'Abd/Add', 'RI/RE'};
    cycle_prefix = 'Cycle #';
    trial = moveck.retrieve_trial(store, output.trial_filenames.task_mocap);
    for i = 1:length(side_names)
        figure('Name', ['Joint kinematics for the ', side_names{i} , ' lower limb']);
        cnt = 1;
        for j = 1:length(lowerlimb_angle_names)
            parameter_name = [side_names{i}, lowerlimb_angle_names{j}];
            prm = trial.retrieve_group(['Analyses/CGM11/Parameters/', parameter_name]);
            prm_experiments = prm.retrieve_set('Experiments').read;
            for k = 1:size(prm_experiments, 1)
                subplot(3,3,cnt);
                hold on;
                num_cycles = size(prm_experiments, 4);
                cycle_names = cell(num_cycles, 1);
                for l = 1:num_cycles
                    plot(0:100, squeeze(prm_experiments(k,1,:,l)));
                    xlim([0,100]);
                    cycle_names{l} = [cycle_prefix, num2str(l)];
                end
                legend(cycle_names);
                ylabel([parameter_name, ' ', angle_suffixes{k}])
                cnt = cnt + 1;
            end
        end
    end
end