function compare_subject_joint_moment(store, output)
    trial = moveck.retrieve_trial(store, output.trial_filenames.task_mocap);
    moment_split_settings = struct();
    moment_split_settings.SourceSet = 'Format/Data/Points';
    moment_split_settings.DestinationGroup = 'References/CGM10/Moments';
    moment_split_settings.Filter = {'Types', 'moment'};
    moment_split_settings.callable_unit = 'data-modifier.set-split';
    if (~trial.exists_group(moment_split_settings.DestinationGroup))
        moveck.transform_data(trial, moment_split_settings);
    end
    references_moment_names = trial.retrieve_group(moment_split_settings.DestinationGroup).list_set_children_name;
    base_moment = {'ShoulderMoment', 'ElbowMoment', 'WristMoment', 'HipMoment', 'KneeMoment', 'AnkleMoment'};
    reference_moment = {};
    computed_moment = {};
    for i = 1:length(base_moment)
        reference_moment{end + 1} = ['L', base_moment{i}];
        reference_moment{end + 1} = ['R', base_moment{i}];
        computed_moment{end + 1} = ['Left', base_moment{i}];
        computed_moment{end + 1} = ['Right', base_moment{i}];
    end
    len_moment_names = length(base_moment) * 2;
    axis_labels = {'X', 'Y', 'Z'};
    fig = figure('NumberTitle','off');
    for a = 1:len_moment_names
        if (~isfield(output.task_analysis, computed_moment{a}))
            continue;
        end
        computed = output.task_analysis.(computed_moment{a});
        reference = trial.retrieve_set([moment_split_settings.DestinationGroup,'/',reference_moment{a}]).read();
        clf(fig)
        set(fig, 'Name', [ '(',num2str(a),'/',num2str(len_moment_names),') ',computed_moment{a}, ' vs ', reference_moment{a} , ' - Space to continue or any key to stop']);
        for i = 1:3
            subplot(3,2,1+2*(i-1));
            hold on;
            data_x = 1:size(computed,1);
            data_yc = computed(:,i);
            data_yr = squeeze(reference(i,1,:)) * output.subject_details.subject_mass;
            plot(data_x, data_yc, 'LineWidth', 2);
            plot(data_x, data_yr, '--', 'LineWidth', 2);  
            title(['Computed (solid) vs reference moment (dashed): ', axis_labels{i}]);
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