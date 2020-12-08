function compare_subject_center_of_mass(store, output)
    trial_mocap = moveck.retrieve_trial(store, output.trial_filenames.task_mocap);
    try
        manufacturer_metadata = trial_mocap.retrieve_group('Format/Metadata/MANUFACTURER');
        if (contains(manufacturer_metadata.retrieve_attribute('SOFTWARE').read(), 'Nexus'))
            version_label = manufacturer_metadata.retrieve_attribute('VERSION_LABEL').read();
            nexus_version = sscanf(version_label, '%i.%i.%s');
            if (nexus_version(1) < 2 || (nexus_version(1) == 2 && nexus_version(2) < 4))
                fprintf('The computation of the Whole body CoM changed in Nexus 2.4 and this trial was generated with Nexus %s\n', nexus_version);
            end
        end
    end
    % Computed set
    model_group = trial_mocap.retrieve_group('Models/CGM11');
    computed_com = model_group.retrieve_set('p_COM').read();
    % Reference set
    reference_com = trial_mocap.retrieve_set('Instruments/Mocap/CentreOfMass').read();
    % Data comparison
    fig = figure('Name', 'Comparison of the Whole Body Center of Mass between CGM-1.0 (dash) and CGM-1.1 (solid)');
    subplot(3,2,1)
    hold on;
    plot(squeeze(computed_com(1,1,:))', 'LineWidth', 2);
    plot(squeeze(reference_com(1,1,:))', '--', 'LineWidth', 2);
    title('Values along the global X axis')
    subplot(3,2,3)
    hold on;
    plot(squeeze(computed_com(2,1,:))', 'LineWidth', 2);
    plot(squeeze(reference_com(2,1,:))', '--', 'LineWidth', 2);
    title('Values along the global Y axis')
    subplot(3,2,5)
    hold on;
    plot(squeeze(computed_com(3,1,:))', 'LineWidth', 2);
    plot(squeeze(reference_com(3,1,:))', '--', 'LineWidth', 2);
    title('Values along the global Z axis')
    subplot(3,2,2)
    plot(squeeze(computed_com(1,1,:) - reference_com(1,1,:))', 'k-', 'LineWidth', 2);
    title('Differences along the X axis')
    subplot(3,2,4)
    plot(squeeze(computed_com(2,1,:) - reference_com(2,1,:))', 'k-', 'LineWidth', 2);
    title('Differences along the Y axis')
    subplot(3,2,6)
    plot(squeeze(computed_com(3,1,:) - reference_com(3,1,:))', 'k-', 'LineWidth', 2);
    title('Differences along the Z axis')
    waitforbuttonpress
    close(fig);
end