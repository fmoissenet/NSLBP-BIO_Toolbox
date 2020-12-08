function plot_segments_pose(store, output, display_options)
    if (nargin == 2)
        display_options = [true, false, false];
    end
    if (length(display_options) < 3)
        display_options(length(display_options)+1:3) = false;
    end
    trial_mocap = store.retrieve_group(output.trial_filenames.task_mocap);
    model_group = trial_mocap.retrieve_group('Models/CGM11');
    figure('Color', [0.1, 0.1, 0.1]);
    hold on;
    axis equal vis3d off
    cameratoolbar
    plot_segment(model_group, 'Head', display_options, [205 92 92] / 256);
    plot_segment(model_group, 'Torso', display_options, [173 255 47] / 256);
    plot_segment(model_group, 'LeftArm', display_options, [255 248 220] / 256);
    plot_segment(model_group, 'LeftForearm', display_options, [255 192 203] / 256);
    plot_segment(model_group, 'LeftHand', display_options, [255 160 122] / 256);
    plot_segment(model_group, 'RightArm', display_options, [255 215 0] / 256);
    plot_segment(model_group, 'RightForearm', display_options, [0 255 255] / 256);
    plot_segment(model_group, 'RightHand', display_options, [255 255 255] / 256);
    plot_segment(model_group, 'Pelvis', display_options, [230 230 250] / 256);
    plot_segment(model_group, 'LeftThigh', display_options, [240 128 128] / 256);
    plot_segment(model_group, 'LeftShank', display_options, [127 255 0] / 256);
    plot_segment(model_group, 'LeftFoot', display_options, [255 235 205] / 256);
    plot_segment(model_group, 'RightThigh', display_options, [255 182 193] / 256);
    plot_segment(model_group, 'RightShank', display_options, [255 127 80] / 256);
    plot_segment(model_group, 'RightFoot', display_options, [255 255 0] / 256);
    view(3);
end

function h = plot_axis(T, axis, scale, opts)
    a = [0, T(1,axis)] * scale + T(1,4);
    b = [0, T(2,axis)] * scale + T(2,4);
    c = [0, T(3,axis)] * scale + T(3,4);
    h = plot3(a, b, c, opts);
end

function h = plot_frame(x,index,scale)
    T = x(:,:,index);
    h(1) = plot_axis(T,1,scale,'r');
    h(2) = plot_axis(T,2,scale,'g');
    h(3) = plot_axis(T,3,scale,'b');
end

function h = plot_segment(group, name, display_options, color)
    % Util functions
    plot3_origin = @(x,opts)plot3(squeeze(x(1,4,:))', squeeze(x(2,4,:))', squeeze(x(3,4,:))', opts{:});
    plot3_com = @(x,opts)plot3(squeeze(x(1,1,:))', squeeze(x(2,1,:))', squeeze(x(3,1,:))', opts{:});
    % - SCS
    pose = group.retrieve_set([name,'/T_SCS']).read();
    % - CoM available?
    is_com_found = false;
    if (group.exists_set([name,'/p_COM']))
        is_com_found = true;
        com = group.retrieve_set([name, '/p_COM']).read();
    end
    num_samples = size(pose, 3);
    sample_indices = 1:10:num_samples;
    if (display_options(1))
        plot3_origin(pose,  {'Color', color, 'LineStyle', '--'})
    end
    if (display_options(2))
        for i = 1:length(sample_indices)
            plot_frame(pose, sample_indices(i), 100);
        end
    end
    if (display_options(3))
        if (is_com_found)
            plot3_com(com, {'o', 'Color', color});
        else
            disp(['The CoM for the segment ' name ' was not found']);
        end
    end
end