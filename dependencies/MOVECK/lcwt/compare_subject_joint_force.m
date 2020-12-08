function compare_subject_joint_force(store, output)
    trial = moveck.retrieve_trial(store, output.trial_filenames.task_mocap);
    force_split_settings = struct();
    force_split_settings.SourceSet = 'Format/Data/Points';
    force_split_settings.DestinationGroup = 'References/CGM10/Forces';
    force_split_settings.Filter = {'Types', 'force'};
    force_split_settings.callable_unit = 'data-modifier.set-split';
    if (~trial.exists_group(force_split_settings.DestinationGroup))
        moveck.transform_data(trial, force_split_settings);
    end
    references_force_names = trial.retrieve_group(force_split_settings.DestinationGroup).list_set_children_name;
    base_force = {'ShoulderForce', 'ElbowForce', 'WristForce', 'HipForce', 'KneeForce', 'AnkleForce'};
    reference_force = {};
    computed_force = {};
    for i = 1:length(base_force)
        reference_force{end + 1} = ['L', base_force{i}];
        reference_force{end + 1} = ['R', base_force{i}];
        computed_force{end + 1} = ['Left', base_force{i}];
        computed_force{end + 1} = ['Right', base_force{i}];
    end
    len_force_names = length(base_force) * 2;
    axis_labels = {'X', 'Y', 'Z'};
    fig = figure('NumberTitle','off');
    for a = 1:len_force_names
        if (~isfield(output.task_analysis, computed_force{a}))
            continue;
        end
        computed = output.task_analysis.(computed_force{a});
        reference = trial.retrieve_set([force_split_settings.DestinationGroup,'/',reference_force{a}]).read();
        clf(fig)
        set(fig, 'Name', [ '(',num2str(a),'/',num2str(len_force_names),') ',computed_force{a}, ' vs ', reference_force{a} , ' - Space to continue or any key to stop']);
        for i = 1:3
            subplot(3,2,1+2*(i-1));
            hold on;
            data_x = 1:size(computed,1);
            data_yc = computed(:,i);
            data_yr = squeeze(reference(i,1,:)) * output.subject_details.subject_mass;
            plot(data_x, data_yc, 'LineWidth', 2);
            plot(data_x, data_yr, '--', 'LineWidth', 2);
            title(['Computed (solid) vs reference force (dashed): ', axis_labels{i}]);
            subplot(3,2,2*i);
            hold on;
            plot(data_x, data_yr - data_yc, 'k-', 'LineWidth', 2)
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