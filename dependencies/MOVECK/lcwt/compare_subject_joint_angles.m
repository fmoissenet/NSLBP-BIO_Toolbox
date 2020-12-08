function compare_subject_joint_angles(store, output)
    trial = moveck.retrieve_trial(store, output.trial_filenames.task_mocap);
    angles_split_settings = struct();
    angles_split_settings.SourceSet = 'Format/Data/Points';
    angles_split_settings.DestinationGroup = 'References/CGM10/Angles';
    angles_split_settings.Filter = {'Types', 'angle'};
    angles_split_settings.callable_unit = 'data-modifier.set-split';
    if (~trial.exists_group(angles_split_settings.DestinationGroup))
        moveck.transform_data(trial, angles_split_settings);
    end
    references_angles_names = trial.retrieve_group(angles_split_settings.DestinationGroup).list_set_children_name;
    base_angles = { ...
        {'HeadAngles', 'HeadProgressionAngles'}, ...
        {'NeckAngles', 'NeckAngles'}, ...
        {'ShoulderAngles', 'ShoulderAngles'}, ...
        {'ElbowAngles', 'ElbowAngles'}, ...
        {'WristAngles', 'WristAngles'}, ...
        {'ThoraxAngles', 'TorsoProgressionAngles'}, ...
        {'SpineAngles', 'SpineAngles'}, ...
        {'PelvisAngles', 'PelvisProgressionAngles'}, ...
        {'HipAngles', 'HipAngles'}, ...
        {'KneeAngles', 'KneeAngles'}, ...
        {'AnkleAngles', 'AnkleAngles'}, ...
        {'FootProgressAngles', 'FootProgressionAngles'}, ...
    };
    reference_angles = {};
    computed_angles = {};
    for i = 1:length(base_angles)
        reference_angles{end + 1} = ['L', base_angles{i}{1}];
        reference_angles{end + 1} = ['R', base_angles{i}{1}];
        computed_angles{end + 1} = ['Left', base_angles{i}{2}];
        computed_angles{end + 1} = ['Right', base_angles{i}{2}];
    end
    len_angles_names = length(base_angles) * 2;
    axis_labels = {'X', 'Y', 'Z'};
    fig = figure('NumberTitle','off');
    for a = 1:len_angles_names
        if (~isfield(output.task_analysis, computed_angles{a}))
            continue;
        end
        computed_angle = output.task_analysis.(computed_angles{a});
        reference_angle = trial.retrieve_set([angles_split_settings.DestinationGroup,'/',reference_angles{a}]).read();
        clf(fig)
        set(fig, 'Name', [ '(',num2str(a),'/',num2str(len_angles_names),') ',computed_angles{a}, ' vs ', reference_angles{a} , ' - Space to continue or any key to stop']);
        for i = 1:3
            subplot(3,2,1+2*(i-1));
            hold on;
            data_x = 1:size(computed_angle,1);
            data_yc = computed_angle(:,i);
            data_yr = squeeze(reference_angle(i,1,:));
            plot(data_x, data_yc, 'LineWidth', 2);
            plot(data_x, data_yr, '--', 'LineWidth', 2);  
            title(['Computed (solid) vs reference angles (dashed): ', axis_labels{i}]);
            subplot(3,2,2*i);
            hold on;
            plot(data_x, data_yr - data_yc, 'k-', 'LineWidth', 2);
            title(['Global difference: ', axis_labels{i}]);
        end
        if (waitforbuttonpress == 1)
            if (get(fig, 'CurrentCharacter') ~= ' ') % Escape character
                break;
            end
        end
    end
    close(fig);
end