function downsample(source, settings)
    settings.callable_unit = 'set-downsample';
    moveck.transform_data(source, settings);
end
