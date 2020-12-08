% Compare (L|R)GroundReaction(Force|Moment) against computed GRW expressed
% at the surface origin
function compare_subject_grw(store, output)
    trial = moveck.retrieve_trial(store, output.trial_filenames.task_mocap);
    reference_group = trial.retrieve_group('References/CGM10');
    processings_group = trial.retrieve_group('Processings');
    if (~processings_group.exists_set('FP1/W_SO'))
        grw_settings = struct();
        grw_settings.callable_unit = 'classical-mechanics.ground-reaction-wrench';
        grw_settings.Location = 'SO';
        grw_settings.SourceGroup = 'Processings/FP1';
        grw_settings.DestinationGroup = 'Processings/FP1';
        moveck.transform_data(trial, grw_settings);
    end
    if (~processings_group.exists_set('FP2/W_SO'))
        grw_settings = struct();
        grw_settings.callable_unit = 'classical-mechanics.ground-reaction-wrench';
        grw_settings.Location = 'SO';
        grw_settings.SourceGroup = 'Processings/FP2';
        grw_settings.DestinationGroup = 'Processings/FP2';
        moveck.transform_data(trial, grw_settings);
    end
    wrenches = struct();
    if (isfield(output.subject_details, 'forceplate_index_left_foot_map'))
        left_wrench = trial.retrieve_set(['Processings/FP',num2str(output.subject_details.forceplate_index_left_foot_map),'/W_SO']).read();
        wrenches.left_force = left_wrench(1:3,:,:);
        wrenches.left_moment = left_wrench(4:6,:,:);
    end
    if (isfield(output.subject_details, 'forceplate_index_right_foot_map'))
        right_wrench = trial.retrieve_set(['Processings/FP',num2str(output.subject_details.forceplate_index_right_foot_map),'/W_SO']).read();
        wrenches.right_force = right_wrench(1:3,:,:);
        wrenches.right_moment = right_wrench(4:6,:,:);
    end
    sets = {...
        {'Forces/LGroundReactionForce',  'left_force'}, ...
        {'Moments/LGroundReactionMoment', 'left_moment'}, ...
        {'Forces/RGroundReactionForce',  'right_force'}, ...
        {'Moments/RGroundReactionMoment', 'right_moment'},...
    };
    len_sets = length(sets);
    axis_labels = {'X', 'Y', 'Z'};
    fig = figure('NumberTitle','off');
    for a = 1:len_sets
        reference = reference_group.retrieve_set(sets{a}{1}).read();
        computed = wrenches.(sets{a}{2});
        clf(fig)
        set(fig, 'Name', [ '(',num2str(a),'/',num2str(len_sets),') ',sets{a}{2}, ' vs ', sets{a}{1} , ' - Space to continue or any key to stop']);
        for i = 1:3
            subplot(3,2,1+2*(i-1));
            hold on;
            inc = floor(size(reference,3) / size(computed,3));
            data_xc = 1:inc:size(reference,3);
            data_yc = squeeze(computed(i,1,:));
            data_xr = 1:size(reference,3);
            data_yr = squeeze(reference(i,1,:)) * output.subject_details.subject_mass;
            plot(data_xc, data_yc, 'LineWidth', 2);
            plot(data_xr, data_yr, '--', 'LineWidth', 2);  
            title(['Computed (solid) vs reference force (dashed): ', axis_labels{i}]);
            subplot(3,2,2*i);
            hold on;
            plot(data_xc, data_yr(data_xc) - data_yc, 'k-', 'LineWidth', 2)
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