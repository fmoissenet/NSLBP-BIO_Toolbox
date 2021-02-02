% Author       : F. Moissenet
%                Kinesiology Laboratory (K-LAB)
%                University of Geneva
%                https://www.unige.ch/medecine/kinesiology
% License      : Creative Commons Attribution-NonCommercial 4.0 International License 
%                https://creativecommons.org/licenses/by-nc/4.0/legalcode
% Source code  : https://github.com/fmoissenet/NSLBP-BIOToolbox
% Reference    : To be defined
% Date         : June 2020
% -------------------------------------------------------------------------
% Description  : To be defined
% Inputs       : To be defined
% Outputs      : To be defined
% -------------------------------------------------------------------------
% Dependencies : - Biomechanical Toolkit (BTK): https://github.com/Biomechanical-ToolKit/BTKCore
% -------------------------------------------------------------------------
% This work is licensed under the Creative Commons Attribution - 
% NonCommercial 4.0 International License. To view a copy of this license, 
% visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% -------------------------------------------------------------------------

function [Trial,tGRF] = ProcessGRFSignals(Session,Trial,GRF,tGRF,fmethod,smethod)

% -------------------------------------------------------------------------
% IDENTIFY FORCEPLATE CYCLES/STEPS
% -------------------------------------------------------------------------

% Right forceplate steps
% -------------------------------------------------------------------------

% Define available stance phase start and stop frames
start = [];
stop  = [];
k1    = 1;
k2    = 1;
for i = 1:size(Trial.Event,2)
    if contains(Trial.Event(i).label,'RHS')
        for j = 1:size(Trial.Event(i).value,2)-1
            start(k1) = Trial.Event(i).value(j);
            k1        = k1+1;
        end
    end
    if contains(Trial.Event(i).label,'RTO')
        for j = 1:size(Trial.Event(i).value,2)
            if ~isempty(start)
                if Trial.Event(i).value(j) > start(1)
                    stop(k2) = Trial.Event(i).value(j);
                    k2       = k2+1;
                end
            end
        end
    end
end

for i = 1:size(start,2)
    
    % Define the foot rectangle during stance
    if ~isempty(Trial.Marker(22).Trajectory.smooth)
        Xmax = max([Trial.Marker(17).Trajectory.smooth(start(i),1) ... % RHEE
                    Trial.Marker(19).Trajectory.smooth(start(i),1) ... % RFMH
                    Trial.Marker(21).Trajectory.smooth(start(i),1) ... % RVMH
                    Trial.Marker(22).Trajectory.smooth(start(i),1) ... % RHAL
                    Trial.Marker(17).Trajectory.smooth(stop(i),1) ...  % RHEE
                    Trial.Marker(19).Trajectory.smooth(stop(i),1) ...  % RFMH
                    Trial.Marker(21).Trajectory.smooth(stop(i),1) ...  % RVMH
                    Trial.Marker(22).Trajectory.smooth(stop(i),1)]);   % RHAL
        Xmin = min([Trial.Marker(17).Trajectory.smooth(start(i),1) ... % RHEE
                    Trial.Marker(19).Trajectory.smooth(start(i),1) ... % RFMH
                    Trial.Marker(21).Trajectory.smooth(start(i),1) ... % RVMH
                    Trial.Marker(22).Trajectory.smooth(start(i),1) ... % RHAL
                    Trial.Marker(17).Trajectory.smooth(stop(i),1) ...  % RHEE
                    Trial.Marker(19).Trajectory.smooth(stop(i),1) ...  % RFMH
                    Trial.Marker(21).Trajectory.smooth(stop(i),1) ...  % RVMH
                    Trial.Marker(22).Trajectory.smooth(stop(i),1)]);   % RHAL
        Ymax = max([Trial.Marker(17).Trajectory.smooth(start(i),2) ... % RHEE
                    Trial.Marker(19).Trajectory.smooth(start(i),2) ... % RFMH
                    Trial.Marker(21).Trajectory.smooth(start(i),2) ... % RVMH
                    Trial.Marker(22).Trajectory.smooth(start(i),2) ... % RHAL
                    Trial.Marker(17).Trajectory.smooth(stop(i),2) ...  % RHEE
                    Trial.Marker(19).Trajectory.smooth(stop(i),2) ...  % RFMH
                    Trial.Marker(21).Trajectory.smooth(stop(i),2) ...  % RVMH
                    Trial.Marker(22).Trajectory.smooth(stop(i),2)]);   % RHAL
        Ymin = min([Trial.Marker(17).Trajectory.smooth(start(i),2) ... % RHEE
                    Trial.Marker(19).Trajectory.smooth(start(i),2) ... % RFMH
                    Trial.Marker(21).Trajectory.smooth(start(i),2) ... % RVMH
                    Trial.Marker(22).Trajectory.smooth(start(i),2) ... % RHAL
                    Trial.Marker(17).Trajectory.smooth(stop(i),2) ...  % RHEE
                    Trial.Marker(19).Trajectory.smooth(stop(i),2) ...  % RFMH
                    Trial.Marker(21).Trajectory.smooth(stop(i),2) ...  % RVMH
                    Trial.Marker(22).Trajectory.smooth(stop(i),2)]);   % RHAL
    else  % Some sessions have no HAL marker, in this case, 20% of SMH-HEE distance is added to the foot
        if mean(Trial.Marker(19).Trajectory.smooth(start(i),1),1) > ... % +X direction
           mean(Trial.Marker(17).Trajectory.smooth(start(i),1),1)
            xmax_offset = 0.2*mean(sqrt((Trial.Marker(17).Trajectory.smooth(:,1)-Trial.Marker(20).Trajectory.smooth(:,1)).^2+...
                                        (Trial.Marker(17).Trajectory.smooth(:,2)-Trial.Marker(20).Trajectory.smooth(:,2)).^2+...
                                        (Trial.Marker(17).Trajectory.smooth(:,3)-Trial.Marker(20).Trajectory.smooth(:,3)).^2),1);
            xmin_offset = Session.markerHeight;
        elseif mean(Trial.Marker(19).Trajectory.smooth(start(i),1),1) < ... % -X direction
               mean(Trial.Marker(17).Trajectory.smooth(start(i),1),1)
            xmax_offset = -Session.markerHeight; 
            xmin_offset = -0.2*mean(sqrt((Trial.Marker(17).Trajectory.smooth(:,1)-Trial.Marker(20).Trajectory.smooth(:,1)).^2+...
                                         (Trial.Marker(17).Trajectory.smooth(:,2)-Trial.Marker(20).Trajectory.smooth(:,2)).^2+...
                                         (Trial.Marker(17).Trajectory.smooth(:,3)-Trial.Marker(20).Trajectory.smooth(:,3)).^2),1);
        end
        Xmax = max([Trial.Marker(17).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(19).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(21).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(17).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(19).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(21).Trajectory.smooth(stop(i),1)]) + xmax_offset;
        Xmin = min([Trial.Marker(17).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(19).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(21).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(17).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(19).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(21).Trajectory.smooth(stop(i),1)]) - xmax_offset;
        Ymax = max([Trial.Marker(17).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(19).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(21).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(17).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(19).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(21).Trajectory.smooth(stop(i),2)]);
        Ymin = min([Trial.Marker(17).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(19).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(21).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(17).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(19).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(21).Trajectory.smooth(stop(i),2)]);
    end
    
    % Check if the foot rectangle is within the forceplate rectangle during
    % the current stance phase
    for j = 1:size(GRF,1)
        if Xmax < Trial.GRF(j).Location.X(2) && ...
           Xmin > Trial.GRF(j).Location.X(1) && ...
           Ymax < Trial.GRF(j).Location.Y(2) && ...
           Ymin > Trial.GRF(j).Location.Y(1)
            Trial.GRF(j).side         = 'R';
            Trial.GRF(j).cycle        = i;
            Trial.GRF(j).Signal.P.raw = GRF(j).P*1e-3; % Convert mm to m
            Trial.GRF(j).Signal.F.raw = GRF(j).F;
            Trial.GRF(j).Signal.M.raw = [zeros(size(GRF(j).M(:,3))) ...
                                         zeros(size(GRF(j).M(:,3))) ...
                                         GRF(j).M(:,3)]*1e-3; % Convert Nmm to Nm
            tGRF(j).P                 = tGRF(j).P*1e-3; % Convert mm to m
            tGRF(j).F                 = tGRF(j).F;
            tGRF(j).M                 = tGRF(j).M*1e-3; % Convert Nmm to Nm
        end
    end
end

% Left forceplate steps
% -------------------------------------------------------------------------

% Define available stance phase start and stop frames
start = [];
stop  = [];
k1    = 1;
k2    = 1;
for i = 1:size(Trial.Event,2)
    if contains(Trial.Event(i).label,'LHS')
        for j = 1:size(Trial.Event(i).value,2)-1
            start(k1) = Trial.Event(i).value(j);
            k1        = k1+1;
        end
    end
    if contains(Trial.Event(i).label,'LTO')
        for j = 1:size(Trial.Event(i).value,2)
            if ~isempty(start)
                if Trial.Event(i).value(j) > start(1)
                    stop(k2) = Trial.Event(i).value(j);
                    k2       = k2+1;
                end
            end
        end
    end
end

for i = 1:size(start,2)
    
    % Define the foot rectangle during stance
    if ~isempty(Trial.Marker(37).Trajectory.smooth)
        Xmax = max([Trial.Marker(32).Trajectory.smooth(start(i),1) ... % RHEE
                    Trial.Marker(34).Trajectory.smooth(start(i),1) ... % RFMH
                    Trial.Marker(36).Trajectory.smooth(start(i),1) ... % RVMH
                    Trial.Marker(37).Trajectory.smooth(start(i),1) ... % RHAL
                    Trial.Marker(32).Trajectory.smooth(stop(i),1) ...  % RHEE
                    Trial.Marker(34).Trajectory.smooth(stop(i),1) ...  % RFMH
                    Trial.Marker(36).Trajectory.smooth(stop(i),1) ...  % RVMH
                    Trial.Marker(37).Trajectory.smooth(stop(i),1)]);   % RHAL
        Xmin = min([Trial.Marker(32).Trajectory.smooth(start(i),1) ... % RHEE
                    Trial.Marker(34).Trajectory.smooth(start(i),1) ... % RFMH
                    Trial.Marker(36).Trajectory.smooth(start(i),1) ... % RVMH
                    Trial.Marker(37).Trajectory.smooth(start(i),1) ... % RHAL
                    Trial.Marker(32).Trajectory.smooth(stop(i),1) ...  % RHEE
                    Trial.Marker(34).Trajectory.smooth(stop(i),1) ...  % RFMH
                    Trial.Marker(36).Trajectory.smooth(stop(i),1) ...  % RVMH
                    Trial.Marker(37).Trajectory.smooth(stop(i),1)]);   % RHAL
        Ymax = max([Trial.Marker(32).Trajectory.smooth(start(i),2) ... % RHEE
                    Trial.Marker(34).Trajectory.smooth(start(i),2) ... % RFMH
                    Trial.Marker(36).Trajectory.smooth(start(i),2) ... % RVMH
                    Trial.Marker(37).Trajectory.smooth(start(i),2) ... % RHAL
                    Trial.Marker(32).Trajectory.smooth(stop(i),2) ...  % RHEE
                    Trial.Marker(34).Trajectory.smooth(stop(i),2) ...  % RFMH
                    Trial.Marker(36).Trajectory.smooth(stop(i),2) ...  % RVMH
                    Trial.Marker(37).Trajectory.smooth(stop(i),2)]);   % RHAL
        Ymin = min([Trial.Marker(32).Trajectory.smooth(start(i),2) ... % RHEE
                    Trial.Marker(34).Trajectory.smooth(start(i),2) ... % RFMH
                    Trial.Marker(36).Trajectory.smooth(start(i),2) ... % RVMH
                    Trial.Marker(37).Trajectory.smooth(start(i),2) ... % RHAL
                    Trial.Marker(32).Trajectory.smooth(stop(i),2) ...  % RHEE
                    Trial.Marker(34).Trajectory.smooth(stop(i),2) ...  % RFMH
                    Trial.Marker(36).Trajectory.smooth(stop(i),2) ...  % RVMH
                    Trial.Marker(37).Trajectory.smooth(stop(i),2)]);   % RHAL
    else % Some sessions have no HAL marker, in this case, 20% of SMH-HEE distance is added to the foot
        if mean(Trial.Marker(34).Trajectory.smooth(start(i),1),1) > ... % +X direction
           mean(Trial.Marker(32).Trajectory.smooth(start(i),1),1)
            xmax_offset = 0.2*mean(sqrt((Trial.Marker(32).Trajectory.smooth(:,1)-Trial.Marker(35).Trajectory.smooth(:,1)).^2+...
                                        (Trial.Marker(32).Trajectory.smooth(:,2)-Trial.Marker(35).Trajectory.smooth(:,2)).^2+...
                                        (Trial.Marker(32).Trajectory.smooth(:,3)-Trial.Marker(35).Trajectory.smooth(:,3)).^2),1);
            xmin_offset = Session.markerHeight;
        elseif mean(Trial.Marker(34).Trajectory.smooth(start(i),1),1) < ... % -X direction
               mean(Trial.Marker(32).Trajectory.smooth(start(i),1),1)
            xmax_offset = -Session.markerHeight; 
            xmin_offset = -0.2*mean(sqrt((Trial.Marker(32).Trajectory.smooth(:,1)-Trial.Marker(35).Trajectory.smooth(:,1)).^2+...
                                         (Trial.Marker(32).Trajectory.smooth(:,2)-Trial.Marker(35).Trajectory.smooth(:,2)).^2+...
                                         (Trial.Marker(32).Trajectory.smooth(:,3)-Trial.Marker(35).Trajectory.smooth(:,3)).^2),1);
        end
        Xmax = max([Trial.Marker(32).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(34).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(36).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(32).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(34).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(36).Trajectory.smooth(stop(i),1)]) + xmax_offset;
        Xmin = min([Trial.Marker(32).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(34).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(36).Trajectory.smooth(start(i),1) ...
                    Trial.Marker(32).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(34).Trajectory.smooth(stop(i),1) ...
                    Trial.Marker(36).Trajectory.smooth(stop(i),1)]) + xmin_offset;
        Ymax = max([Trial.Marker(32).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(34).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(36).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(32).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(34).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(36).Trajectory.smooth(stop(i),2)]);
        Ymin = min([Trial.Marker(32).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(34).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(36).Trajectory.smooth(start(i),2) ...
                    Trial.Marker(32).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(34).Trajectory.smooth(stop(i),2) ...
                    Trial.Marker(36).Trajectory.smooth(stop(i),2)]);
    end
    
    % Check if the foot rectangle is within the forceplate rectangle during
    % the current stance phase
    for j = 1:size(GRF,1)
        if Xmax < Trial.GRF(j).Location.X(2) && ...
           Xmin > Trial.GRF(j).Location.X(1) && ...
           Ymax < Trial.GRF(j).Location.Y(2) && ...
           Ymin > Trial.GRF(j).Location.Y(1)
            Trial.GRF(j).side         = 'L';
            Trial.GRF(j).cycle        = i;
            Trial.GRF(j).Signal.P.raw = GRF(j).P*1e-3; % Convert mm to m
            Trial.GRF(j).Signal.F.raw = GRF(j).F;
            Trial.GRF(j).Signal.M.raw = [zeros(size(GRF(j).M(:,3))) ...
                                         zeros(size(GRF(j).M(:,3))) ...
                                         GRF(j).M(:,3)]*1e-3; % Convert Nmm to Nm
            tGRF(j).P                 = tGRF(j).P*1e-3; % Convert mm to m
            tGRF(j).F                 = tGRF(j).F;
            tGRF(j).M                 = tGRF(j).M*1e-3; % Convert Nmm to Nm
        end
    end
end

% -------------------------------------------------------------------------
% SIGNAL FILTERING
% -------------------------------------------------------------------------
if ~isempty(Trial.GRF)
    for i = 1:size(Trial.GRF,2)
        
        % Method 1: No filtering
        if strcmp(fmethod.type,'none') 
            Trial.GRF(i).Signal.P.filt   = Trial.GRF(i).Signal.P.raw;
            Trial.GRF(i).Signal.F.filt   = Trial.GRF(i).Signal.F.raw;
            Trial.GRF(i).Signal.M.filt   = Trial.GRF(i).Signal.M.raw;
            Trial.GRF(i).Processing.filt = fmethod.type;
        
        % Method 2: Vertical force threshold ([fmethod.parameter] N)
        elseif strcmp(fmethod.type,'threshold') 
            if ~isempty(Trial.GRF(i).Signal.F.raw)
                for j = 1:size(Trial.GRF(i).Signal.F.raw,1)
                    if Trial.GRF(i).Signal.F.raw(j,3) < fmethod.parameter
                        Trial.GRF(i).Signal.P.filt(j,:) = zeros(1,3);
                        Trial.GRF(i).Signal.F.filt(j,:) = zeros(1,3);
                        Trial.GRF(i).Signal.M.filt(j,:) = zeros(1,3);
                        tGRF(i).P(j,:) = zeros(1,3);
                        tGRF(i).F(j,:) = zeros(1,3);
                        tGRF(i).M(j,:) = zeros(1,3);
                    else
                        Trial.GRF(i).Signal.P.filt(j,:) = Trial.GRF(i).Signal.P.raw(j,:);
                        Trial.GRF(i).Signal.F.filt(j,:) = Trial.GRF(i).Signal.F.raw(j,:);
                        Trial.GRF(i).Signal.M.filt(j,:) = Trial.GRF(i).Signal.M.raw(j,:);
                        tGRF(i).P(j,:) = tGRF(i).P(j,:);
                        tGRF(i).F(j,:) = tGRF(i).F(j,:);
                        tGRF(i).M(j,:) = tGRF(i).M(j,:);
                    end
                end
            end
            Trial.GRF(i).Processing.filt = fmethod.type;
        end
        
    end
end

% -------------------------------------------------------------------------
% SIGNAL SMOOTHING
% -------------------------------------------------------------------------
if ~isempty(Trial.GRF)
    for i = 1:size(Trial.GRF,2)

        % Method 1: No smoothing
        if strcmp(smethod.type,'none') 
            if ~isempty(Trial.GRF(i).Signal.F.raw)
                Trial.GRF(i).Signal.P.smooth = Trial.GRF(i).Signal.P.filt;
                Trial.GRF(i).Signal.F.smooth = Trial.GRF(i).Signal.F.filt;
                Trial.GRF(i).Signal.M.smooth = Trial.GRF(i).Signal.M.filt;
            end
            Trial.GRF(i).Processing.smooth = smethod.type;

        % Method 2: Low pass filter (Butterworth 2nd order, [smethod.parameter] Hz)
        elseif strcmp(smethod.type,'butterLow2') 
            if ~isempty(Trial.GRF(i).Signal.F.raw)
                [B,A]                          = butter(1,smethod.parameter/(Trial.fanalog/2),'low');
                Trial.GRF(i).Signal.P.smooth   = filtfilt(B,A,Trial.GRF(i).Signal.P.filt);
                Trial.GRF(i).Signal.F.smooth   = filtfilt(B,A,Trial.GRF(i).Signal.F.filt);
                Trial.GRF(i).Signal.M.smooth   = filtfilt(B,A,Trial.GRF(i).Signal.M.filt);
                tGRF(i).P                      = filtfilt(B,A,tGRF(i).P);
                tGRF(i).F                      = filtfilt(B,A,tGRF(i).F);
                tGRF(i).M                      = filtfilt(B,A,tGRF(i).M);
            end
            Trial.GRF(i).Processing.smooth = smethod.type;
        end

    end
end