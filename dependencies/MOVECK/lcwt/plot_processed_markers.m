function plot_processed_markers(trial)
    mocap_instrument = trial.retrieve_group('Instruments/Mocap');
    mocap_processing = trial.retrieve_group('Processings/Mocap');
    marker_names = mocap_instrument.list_set_children_name();
    len_marker_names = length(marker_names);
    axis_labels = {'X', 'Y', 'Z'};
    fig = figure('NumberTitle','off');
    for i = 1:len_marker_names
        markers_raw = mocap_instrument.retrieve_set(marker_names{i}).read();
        markers_processed = mocap_processing.retrieve_set(marker_names{i}).read();
        clf(fig)
        set(fig, 'Name', [ '(',num2str(i),'/',num2str(len_marker_names),') ',marker_names{i}, ' - Space to continue or any key to stop']);
        for i = 1:3
            subplot(3,2,1+2*(i-1));
            hold on;
            data_x = 1:size(markers_raw,3);
            data_yr = squeeze(markers_raw(i,1,:))';
            data_yp = squeeze(markers_processed(i,1,:))';
            plot(data_x, data_yr);
            nan_raw = isnan(data_yr);
            plot(data_x(nan_raw), repmat(min(data_yr), length(data_x(nan_raw)), 1), 'ks');
            plot(data_x, data_yp,'--');
            nan_processed = isnan(data_yp);
            plot(data_x(nan_processed), repmat(min(data_yp), length(data_x(nan_processed)), 1), 'kx');
            title(['Raw vs Processed data: ', axis_labels{i}]);
            subplot(3,2,2*i);
            hold on;
            plot(data_x, data_yp - data_yr)
            plot(data_x(nan_raw), zeros(length(data_x(nan_raw)), 1), 'ks');
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
