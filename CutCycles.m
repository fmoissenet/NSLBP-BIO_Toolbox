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
    end 
    
else % Movements other than gait
            
    if contains(Trial.Event.label,'start')

        % Set temporal parameters
        kr = 1;
        for j = 1:size(Trial.Event.value,2)-1
            % Markers
            Cycle(kr).start  = Trial.Event.value(j);
            Cycle(kr).stop   = Trial.Event.value(j+1);
            Cycle(kr).n      = Cycle(kr).stop-Cycle(kr).start+1;
            Cycle(kr).k      = (1:Cycle(kr).n)';
            Cycle(kr).k0     = (linspace(1,Cycle(kr).n,101))';
            % Analogs
            Cycle(kr).starta = Trial.Event.value(j)*Trial.fanalog/Trial.fmarker;
            Cycle(kr).stopa  = Trial.Event.value(j+1)*Trial.fanalog/Trial.fmarker;
            Cycle(kr).na     = Cycle(kr).stopa-Cycle(kr).starta+1;
            Cycle(kr).ka     = (1:Cycle(kr).na)';
            Cycle(kr).k0a    = (linspace(1,Cycle(kr).na,101))';
            kr                = kr+1;
        end               
    end
    
    % Cut cycles
    for i = 1:size(Cycle,2)
        
        % Marker trajectories
        for j = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(j).Trajectory.smooth)
                temp = interp1(Cycle(i).k,...
                               Trial.Marker(j).Trajectory.smooth(Cycle(i).start:Cycle(i).stop,:),...
                               Cycle(i).k0,...
                               'spline');
                Trial.Marker(j).Trajectory.cycle(:,:,i) = temp;
                clear temp;
            end
        end
        
        % EMG signals
        for j = 1:size(Trial.EMG,2)
            if ~isempty(Trial.EMG(j).Signal.smooth)
                temp = interp1(Cycle(i).ka,...
                               Trial.EMG(j).Signal.smooth(Cycle(i).starta:Cycle(i).stopa,:),...
                               Cycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.cycle(:,:,i) = temp;
                clear temp;
            end
            if ~isempty(Trial.EMG(j).Signal.norm)
                temp = interp1(Cycle(i).ka,...
                               Trial.EMG(j).Signal.norm(Cycle(i).starta:Cycle(i).stopa,:),...
                               Cycle(i).k0a,...
                               'spline');
                Trial.EMG(j).Signal.cyclen(:,:,i) = temp;
                clear temp;
            end
        end
        
        % GRF signals
        for j = 1:size(Trial.GRF,2)
            if ~isempty(Trial.GRF(j).Signal.P.smooth)
                temp = interp1(Cycle(i).ka,...
                               Trial.GRF(j).Signal.P.smooth(Cycle(i).starta:Cycle(i).stopa,:),...
                               Cycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.P.rcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(Cycle(i).ka,...
                               Trial.GRF(j).Signal.F.smooth(Cycle(i).starta:Cycle(i).stopa,:),...
                               Cycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.F.rcycle(:,:,i) = temp;
                clear temp;
                temp = interp1(Cycle(i).ka,...
                               Trial.GRF(j).Signal.M.smooth(Cycle(i).starta:Cycle(i).stopa,:),...
                               Cycle(i).k0a,...
                               'spline');
                Trial.GRF(j).Signal.M.rcycle(:,:,i) = temp;
                clear temp;
            end
        end
    end
end