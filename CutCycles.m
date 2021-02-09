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
% Description  : This routine aims to process 3D marker trajectories
% Inputs       : To be defined
% Outputs      : To be defined
% -------------------------------------------------------------------------
% Dependencies : To be defined
% -------------------------------------------------------------------------
% This work is licensed under the Creative Commons Attribution - 
% NonCommercial 4.0 International License. To view a copy of this license, 
% visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% -------------------------------------------------------------------------

function Trial = CutCycles(Trial)

RCycle = [];
LCycle = [];
Cycle = [];

% Walking trials
if contains(Trial.type,'Gait')
    
    % Identify cycles
    for i = 1:size(Trial.Event,2)
        
        % Right gait cycles
        if contains(Trial.Event(i).label,'RHS')
        
            % Set temporal parameters
            kr = 1;
            for j = 1:size(Trial.Event(i).value,2)-1
                % Markers
                RCycle(kr).start  = Trial.Event(i).value(j);
                RCycle(kr).stop   = Trial.Event(i).value(j+1);
                RCycle(kr).n      = RCycle(kr).stop-RCycle(kr).start+1;
                RCycle(kr).k      = (1:RCycle(kr).n)';
                RCycle(kr).k0     = (linspace(1,RCycle(kr).n,101))';
                % Analogs
                RCycle(kr).starta = Trial.Event(i).value(j)*Trial.fanalog/Trial.fmarker;
                RCycle(kr).stopa  = Trial.Event(i).value(j+1)*Trial.fanalog/Trial.fmarker;
                RCycle(kr).na     = RCycle(kr).stopa-RCycle(kr).starta+1;
                RCycle(kr).ka     = (1:RCycle(kr).na)';
                RCycle(kr).k0a    = (linspace(1,RCycle(kr).na,101))';
                kr                = kr+1;
            end

        % Left gait cycles
        elseif contains(Trial.Event(i).label,'LHS')
        
            % Set temporal parameters
            kl = 1;
            for j = 1:size(Trial.Event(i).value,2)-1
                % Markers
                LCycle(kl).start  = Trial.Event(i).value(j);
                LCycle(kl).stop   = Trial.Event(i).value(j+1);
                LCycle(kl).n      = LCycle(kl).stop-LCycle(kl).start+1;
                LCycle(kl).k      = (1:LCycle(kl).n)';
                LCycle(kl).k0     = (linspace(1,LCycle(kl).n,101))';
                % Analogs                
                LCycle(kl).starta = Trial.Event(i).value(j)*Trial.fanalog/Trial.fmarker;
                LCycle(kl).stopa  = Trial.Event(i).value(j+1)*Trial.fanalog/Trial.fmarker;
                LCycle(kl).na     = LCycle(kl).stopa-LCycle(kl).starta+1;
                LCycle(kl).ka     = (1:LCycle(kl).na)';
                LCycle(kl).k0a    = (linspace(1,LCycle(kl).na,101))';
                kl                = kl+1;
            end
        end
    end
    
    % Cut right gait cycles
    for i = 1:size(RCycle,2)
        
        % Marker trajectories
        for j = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(j).Trajectory.smooth)
                temp = interp1(RCycle(i).k,...
                               Trial.Marker(j).Trajectory.smooth(RCycle(i).start:RCycle(i).stop,:),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Marker(j).Trajectory.rcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % Vmarker trajectories
        for j = 1:size(Trial.Vmarker,2)
            if ~isempty(Trial.Vmarker(j).Trajectory.smooth)
                temp = interp1(RCycle(i).k,...
                               Trial.Vmarker(j).Trajectory.smooth(RCycle(i).start:RCycle(i).stop,:),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Vmarker(j).Trajectory.rcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % EMG signals
        for j = 1:size(Trial.EMG,2)
            if ~isempty(Trial.EMG(j).Signal.smooth)
                temp = interp1(RCycle(i).ka,...
                               Trial.EMG(j).Signal.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.rcycle(:,:,i) = temp;
                clear temp;
            end
            if ~isempty(Trial.EMG(j).Signal.norm)
                temp = interp1(RCycle(i).ka,...
                               Trial.EMG(j).Signal.norm(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.rcyclen(:,:,i) = temp;
                clear temp;
            end
        end
        
        % GRF signals
        for j = 1:size(Trial.GRF,2)
            if ~isempty(Trial.GRF(j).Signal.P.smooth)
                temp = interp1(RCycle(i).ka,...
                               Trial.GRF(j).Signal.P.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.P.rcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).ka,...
                               Trial.GRF(j).Signal.F.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.F.rcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).ka,...
                               Trial.GRF(j).Signal.M.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.M.rcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % Segment kinematics
        for j = 1:size(Trial.Segment,2)
            if ~isempty(Trial.Segment(j).T.smooth)
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).rM.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).rM.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).Q.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Q.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).T.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).T.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).Euler.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Euler.rcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
        
        % Joint kinematics
        for j = 1:size(Trial.Joint,2)
            if ~isempty(Trial.Joint(j).T.smooth)
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Joint(j).T.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Joint(j).T.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Joint(j).Euler.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Joint(j).Euler.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Joint(j).dj.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Joint(j).dj.rcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
    end
    
    % Cut left gait cycles
    for i = 1:size(LCycle,2)
        
        % Marker trajectories
        for j = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(j).Trajectory.smooth)
                temp = interp1(LCycle(i).k,...
                               Trial.Marker(j).Trajectory.smooth(LCycle(i).start:LCycle(i).stop,:),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Marker(j).Trajectory.lcycle(:,:,i) = temp;
                clear temp;
                
            end
        end
        
        % Vmarker trajectories
        for j = 1:size(Trial.Vmarker,2)
            if ~isempty(Trial.Vmarker(j).Trajectory.smooth)
                temp = interp1(LCycle(i).k,...
                               Trial.Vmarker(j).Trajectory.smooth(LCycle(i).start:LCycle(i).stop,:),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Vmarker(j).Trajectory.lcycle(:,:,i) = temp;
                clear temp;
                
            end
        end
        
        % EMG signals
        for j = 1:size(Trial.EMG,2)
            if ~isempty(Trial.EMG(j).Signal.smooth)
                temp = interp1(LCycle(i).ka,...
                               Trial.EMG(j).Signal.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.lcycle(:,:,i) = temp;
                clear temp;
            end
            if ~isempty(Trial.EMG(j).Signal.norm)
                temp = interp1(LCycle(i).ka,...
                               Trial.EMG(j).Signal.norm(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.lcyclen(:,:,i) = temp;
                clear temp;
            end
        end
        
        % GRF signals
        for j = 1:size(Trial.GRF,2)
            if ~isempty(Trial.GRF(j).Signal.P.smooth)
                temp = interp1(LCycle(i).ka,...
                               Trial.GRF(j).Signal.P.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.P.lcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).ka,...
                               Trial.GRF(j).Signal.F.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.F.lcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).ka,...
                               Trial.GRF(j).Signal.M.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.M.lcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % Segment kinematics
        for j = 1:size(Trial.Segment,2)
            if ~isempty(Trial.Segment(j).T.smooth)
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).rM.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).rM.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).Q.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Q.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).T.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).T.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).Euler.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Euler.lcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
        
        % Joint kinematics
        for j = 1:size(Trial.Joint,2)
            if ~isempty(Trial.Joint(j).T.smooth)
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Joint(j).T.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Joint(j).T.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Joint(j).Euler.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Joint(j).Euler.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Joint(j).dj.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Joint(j).dj.lcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
    end 

% Other movements
elseif contains(Trial.type,'S2S') || ...
       contains(Trial.type,'Trunk') || ...
       contains(Trial.type,'Weight') || ...
       contains(Trial.type,'Perturbation')

    % Identify rcycle, i.e. Part 1 of the full cycle (start to back)
    for i = 1:size(Trial.Event,2)
        if contains(Trial.Event(i).label,'start')
            % Set temporal parameters
            kr = 1;
            for j = 1:size(Trial.Event(i).value,2)-1
                % Markers
                RCycle(kr).start  = Trial.Event(i).value(j);
                % Analogs
                RCycle(kr).starta = Trial.Event(i).value(j)*Trial.fanalog/Trial.fmarker;
                kr                = kr+1;
            end       
        elseif contains(Trial.Event(i).label,'back')
            % Set temporal parameters
            kr = 1;
            for j = 1:size(Trial.Event(i).value,2)
                if Trial.Event(i).value(j) > RCycle(kr).start
                    % Markers
                    RCycle(kr).stop  = Trial.Event(i).value(j);
                    RCycle(kr).n     = RCycle(kr).stop-RCycle(kr).start+1;
                    RCycle(kr).k     = (1:RCycle(kr).n)';
                    RCycle(kr).k0    = (linspace(1,RCycle(kr).n,101))';
                    % Analogs
                    RCycle(kr).stopa = Trial.Event(i).value(j)*Trial.fanalog/Trial.fmarker;
                    RCycle(kr).na    = RCycle(kr).stopa-RCycle(kr).starta+1;
                    RCycle(kr).ka    = (1:RCycle(kr).na)';
                    RCycle(kr).k0a   = (linspace(1,RCycle(kr).na,101))';
                    kr               = kr+1;
                end
            end
        end
    end
    
    % Identify lcycle, i.e. Part 2 of the full cycle (start to back)
    for i = 1:size(Trial.Event,2)
        if contains(Trial.Event(i).label,'start')
            % Set temporal parameters
            kl = 1;
            for j = 2:size(Trial.Event(i).value,2)
                % Markers
                LCycle(kl).stop  = Trial.Event(i).value(j);
                % Analogs
                LCycle(kl).stopa = Trial.Event(i).value(j)*Trial.fanalog/Trial.fmarker;
                kl               = kl+1;
            end 
        elseif contains(Trial.Event(i).label,'back')
            % Set temporal parameters
            kl = 1;
            for j = 1:size(Trial.Event(i).value,2)
                % Markers
                LCycle(kl).start = Trial.Event(i).value(j);
                LCycle(kl).n     = LCycle(kl).stop-LCycle(kl).start+1;
                LCycle(kl).k     = (1:LCycle(kl).n)';
                LCycle(kl).k0    = (linspace(1,LCycle(kl).n,101))';
                % Analogs
                LCycle(kl).starta = Trial.Event(i).value(j)*Trial.fanalog/Trial.fmarker;
                LCycle(kl).na     = LCycle(kl).stopa-LCycle(kl).starta+1;
                LCycle(kl).ka     = (1:LCycle(kl).na)';
                LCycle(kl).k0a    = (linspace(1,LCycle(kl).na,101))';
                kl                = kl+1;
            end  
        end
    end
    
    % Cut right gait cycles
    for i = 1:size(RCycle,2)
        
        % Marker trajectories
        for j = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(j).Trajectory.smooth)
                temp = interp1(RCycle(i).k,...
                               Trial.Marker(j).Trajectory.smooth(RCycle(i).start:RCycle(i).stop,:),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Marker(j).Trajectory.rcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % Vmarker trajectories
        for j = 1:size(Trial.Vmarker,2)
            if ~isempty(Trial.Vmarker(j).Trajectory.smooth)
                temp = interp1(RCycle(i).k,...
                               Trial.Vmarker(j).Trajectory.smooth(RCycle(i).start:RCycle(i).stop,:),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Vmarker(j).Trajectory.rcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % EMG signals
        for j = 1:size(Trial.EMG,2)
            if ~isempty(Trial.EMG(j).Signal.smooth)
                temp = interp1(RCycle(i).ka,...
                               Trial.EMG(j).Signal.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.rcycle(:,:,i) = temp;
                clear temp;
            end
            if ~isempty(Trial.EMG(j).Signal.norm)
                temp = interp1(RCycle(i).ka,...
                               Trial.EMG(j).Signal.norm(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.rcyclen(:,:,i) = temp;
                clear temp;
            end
        end
        
        % GRF signals
        for j = 1:size(Trial.GRF,2)
            if ~isempty(Trial.GRF(j).Signal.P.smooth)
                temp = interp1(RCycle(i).ka,...
                               Trial.GRF(j).Signal.P.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.P.rcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).ka,...
                               Trial.GRF(j).Signal.F.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.F.rcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).ka,...
                               Trial.GRF(j).Signal.M.smooth(RCycle(i).starta:RCycle(i).stopa,:),...
                               RCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.M.rcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % Segment kinematics
        for j = 1:size(Trial.Segment,2)
            if ~isempty(Trial.Segment(j).T.smooth)
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).rM.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).rM.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).Q.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Q.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).T.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).T.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Segment(j).Euler.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Euler.rcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
        
        % Joint kinematics
        for j = 1:size(Trial.Joint,2)
            if ~isempty(Trial.Joint(j).T.smooth)
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Joint(j).T.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Joint(j).T.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Joint(j).Euler.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Joint(j).Euler.rcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(RCycle(i).k,...
                               permute(Trial.Joint(j).dj.smooth(:,:,RCycle(i).start:RCycle(i).stop),[3,1,2]),...
                               RCycle(i).k0,...
                               'spline');
                Trial.Joint(j).dj.rcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
    end
    
    % Cut left gait cycles
    for i = 1:size(LCycle,2)
        
        % Marker trajectories
        for j = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(j).Trajectory.smooth)
                temp = interp1(LCycle(i).k,...
                               Trial.Marker(j).Trajectory.smooth(LCycle(i).start:LCycle(i).stop,:),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Marker(j).Trajectory.lcycle(:,:,i) = temp;
                clear temp;
                
            end
        end
        
        % Vmarker trajectories
        for j = 1:size(Trial.Vmarker,2)
            if ~isempty(Trial.Vmarker(j).Trajectory.smooth)
                temp = interp1(LCycle(i).k,...
                               Trial.Vmarker(j).Trajectory.smooth(LCycle(i).start:LCycle(i).stop,:),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Vmarker(j).Trajectory.lcycle(:,:,i) = temp;
                clear temp;
                
            end
        end
        
        % EMG signals
        for j = 1:size(Trial.EMG,2)
            if ~isempty(Trial.EMG(j).Signal.smooth)
                temp = interp1(LCycle(i).ka,...
                               Trial.EMG(j).Signal.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.lcycle(:,:,i) = temp;
                clear temp;
            end
            if ~isempty(Trial.EMG(j).Signal.norm)
                temp = interp1(LCycle(i).ka,...
                               Trial.EMG(j).Signal.norm(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.lcyclen(:,:,i) = temp;
                clear temp;
            end
        end
        
        % GRF signals
        for j = 1:size(Trial.GRF,2)
            if ~isempty(Trial.GRF(j).Signal.P.smooth)
                temp = interp1(LCycle(i).ka,...
                               Trial.GRF(j).Signal.P.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.P.lcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).ka,...
                               Trial.GRF(j).Signal.F.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.F.lcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).ka,...
                               Trial.GRF(j).Signal.M.smooth(LCycle(i).starta:LCycle(i).stopa,:),...
                               LCycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.M.lcycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % Segment kinematics
        for j = 1:size(Trial.Segment,2)
            if ~isempty(Trial.Segment(j).T.smooth)
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).rM.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).rM.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).Q.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Q.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).T.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).T.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Segment(j).Euler.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Segment(j).Euler.lcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
        
        % Joint kinematics
        for j = 1:size(Trial.Joint,2)
            if ~isempty(Trial.Joint(j).T.smooth)
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Joint(j).T.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Joint(j).T.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Joint(j).Euler.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Joint(j).Euler.lcycle(:,:,:,i) = temp;
                clear temp;
                temp = interp1(LCycle(i).k,...
                               permute(Trial.Joint(j).dj.smooth(:,:,LCycle(i).start:LCycle(i).stop),[3,1,2]),...
                               LCycle(i).k0,...
                               'spline');
                Trial.Joint(j).dj.lcycle(:,:,:,i) = temp;
                clear temp;
            end
        end
    end 
end