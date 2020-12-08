function create_instrument(source, settings)
    settings.callable_unit = 'instrument-build';
    moveck.transform_data(source, settings);
end