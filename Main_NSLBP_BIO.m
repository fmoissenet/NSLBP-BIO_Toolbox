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
% Description  : Main routine used to launch NSLBP-BIO routines
% -------------------------------------------------------------------------
% Dependencies : To be defined
% -------------------------------------------------------------------------
% This work is licensed under the Creative Commons Attribution - 
% NonCommercial 4.0 International License. To view a copy of this license, 
% visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% INIT THE WORKSPACE
% -------------------------------------------------------------------------
clearvars;
close all;
clc;

% -------------------------------------------------------------------------
% SET FOLDERS
% -------------------------------------------------------------------------
disp('Set folders');
Folder.toolbox      = 'C:\Users\moissene\Documents\Professionnel\projets recherche\2019 - NSCLBP - Biomarkers\Données\NSLBP-BIO_Toolbox\';
Folder.dependencies = [Folder.toolbox,'dependencies\'];
Folder.data         = uigetdir;
Folder.export       = [Folder.data,'\output\'];
Folder.biomarkers   = 'C:\Users\moissene\Documents\Professionnel\projets recherche\2019 - NSCLBP - Biomarkers\Données\';
addpath(Folder.toolbox);
addpath(genpath(Folder.dependencies));

% -------------------------------------------------------------------------
% DEFINE PARTICIPANT
% -------------------------------------------------------------------------
disp('Set participant parameters');
Participant.id           = '001';
Participant.type         = 'Control';
Participant.gender       = 'Male';
Participant.inclusionAge = 27; % years
Participant.pelvisWidth  = 27; % cm
Participant.RLegLength   = 93; % cm
Participant.LLegLength   = 93; % cm
Participant.RKneeWidth   = 11; % cm
Participant.LKneeWidth   = 11; % cm
Participant.RAnkleWidth  = 7; % cm
Participant.LAnkleWidth  = 7; % cm

% -------------------------------------------------------------------------
% DEFINE SESSION
% -------------------------------------------------------------------------
disp('Set session parameters');
Session.date              = '03-06-2020';
Session.type              = 'INI_session';
Session.examiner          = 'FM';
Session.participantHeight = 180.5; % cm
Session.participantWeight = 81.0; % kg
Session.markerHeight      = 14*1e-3; % m

% -------------------------------------------------------------------------
% LOAD C3D FILES
% -------------------------------------------------------------------------
disp('Extract data from C3D files');

% List all trial types
trialTypes = {'Static',...
              'Endurance_Ito','Endurance_Sorensen',...
              'Gait_Fast','Gait_Normal','Gait_Slow',...
              'Posture_Standing','Posture_Sitting',...
              'Perturbation_R_Shoulder','Perturbation_L_Shoulder',...
              'S2S_Constrained','S2S_Unconstrained',...
              'Swing_R_Leg','Swing_L_Leg',...
              'Trunk_Forward','Trunk_Lateral','Trunk_Rotation',...
              'Weight_Constrained','Weight_Unconstrained',...
              'sMVC'};

% Extract data from C3D files
cd(Folder.data);
c3dFiles = dir('*.c3d');
k1       = 1;
k2       = 1;
for i = 1:size(c3dFiles)
    disp(['  - ',c3dFiles(i).name]);
    for j = 1:size(trialTypes,2)
        if strfind(c3dFiles(i).name,trialTypes{j}) == 1
            if ~isempty(strfind(trialTypes{j},'Static'))
                Static(k1).type    = trialTypes{j};
                Static(k1).file    = c3dFiles(i).name;
                Static(k1).btk     = btkReadAcquisition(c3dFiles(i).name);
                Static(k1).n0      = btkGetFirstFrame(Static(k1).btk);
                Static(k1).n1      = btkGetLastFrame(Static(k1).btk)-Static(k1).n0+1;
                Static(k1).fmarker = btkGetPointFrequency(Static(k1).btk);
                Static(k1).fanalog = btkGetAnalogFrequency(Static(k1).btk);
                k1 = k1+1;
            else
                Trial(k2).type    = trialTypes{j};
                Trial(k2).file    = c3dFiles(i).name;
                Trial(k2).btk     = btkReadAcquisition(c3dFiles(i).name);
                Trial(k2).n0      = btkGetFirstFrame(Trial(k2).btk);
                Trial(k2).n1      = btkGetLastFrame(Trial(k2).btk)-Trial(k2).n0+1;
                Trial(k2).fmarker = btkGetPointFrequency(Trial(k2).btk);
                Trial(k2).fanalog = btkGetAnalogFrequency(Trial(k2).btk);
                k2 = k2+1;
            end
        end
    end
end
clear k1 k2 c3dFiles trialTypes;

% -------------------------------------------------------------------------
% PRE-PROCESS DATA
% -------------------------------------------------------------------------

% Static data
disp('Pre-process static data');
for i = 1:size(Static,2)
    disp(['  - ',Static(i).file]);
    
    % Get manually defined events
    Static(i).Event = [];
    
    % Process marker trajectories
    Marker           = btkGetMarkers(Static(i).btk);
    Static(i).Marker = [];
    Static(i)        = InitialiseMarkerTrajectories(Static(i),Marker);
    Static(i)        = ProcessMarkerTrajectories([],Static(i));
    clear Marker;
    
    % Process EMG signals
    Static(i).EMG = [];
    
    % Process forceplate signals
    Static(i).GRF = [];
    
    % Store processed static data in a new C3D file
    ExportC3D(Static(i),[],Participant,Session,Folder);
end

% EMG calibration data
disp('Pre-process EMG calibration data');
Calibration(1).type = 'EMG_calibration';
Calibration(1).EMG  = [];
for i = 1:size(Trial,2)
    
    if ~isempty(strfind(Trial(i).type,'Endurance')) || ...
       ~isempty(strfind(Trial(i).type,'sMVC')) 
        
        disp(['  - ',Trial(i).file]);

        % Get manually defined events
        Event          = btkGetEvents(Trial(i).btk);
        Trial(i).Event = [];
        Trial(i)       = InitialiseEvents(Trial(i),Event);
        clear Event;   
    
        % Process marker trajectories
        Trial(i).Marker = [];  
        
        % Process EMG signals
        EMG                    = btkGetAnalogs(Trial(i).btk);
        Trial(i).EMG           = [];
        Trial(i)               = InitialiseEMGSignals(Trial(i),EMG);
        fmethod.type           = 'butterBand4';
        fmethod.parameter      = [10 450];
        smethod.type           = 'butterLow2';
        smethod.parameter      = 3;
        nmethod.type           = 'sMVC';
        [Calibration,Trial(i)] = ProcessEMGSignals(Calibration,Trial(i),1,fmethod,smethod,nmethod);
        clear EMG fmethod smethod nmethod;
    
        % Process forceplate signals
        Trial(i).GRF = [];
        
        % Store processed static data in a new C3D file
        ExportC3D(Trial(i),[],Participant,Session,Folder);        
   
    end
end

% Trial data
disp('Pre-process trial data');
for i = 1:size(Trial,2)
    
    if isempty(strfind(Trial(i).type,'sMVC')) % Endurance tasks considered as Trial here

        disp(['  - ',Trial(i).file]);

        % Get manually defined events
        Trial(i).Event = [];
        Event          = btkGetEvents(Trial(i).btk);
        Trial(i)       = InitialiseEvents(Trial(i),Event);
        clear Event;   

        % Process marker trajectories   
        Trial(i).Marker      = [];
        Marker               = btkGetMarkers(Trial(i).btk);
        Trial(i)             = InitialiseMarkerTrajectories(Trial(i),Marker);        
        fmethod.type         = 'intercor';
        fmethod.gapThreshold = [];
        smethod.type         = 'movmean';
        smethod.parameter    = 15;        
        Trial(i)             = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);        
        clear Marker fmethod smethod;
        
        % Compute segment and joint kinematics
        Trial(i).Vmarker = [];
        Trial(i).Segment = [];
        Trial(i).Joint   = [];
        Trial(i)         = InitialiseVmarkerTrajectories(Trial(i));
        Trial(i)         = InitialiseSegments(Trial(i));
        Trial(i)         = InitialiseJoints(Trial(i));
        if isempty(strfind(Trial(i).type,'Endurance'))
            Trial(i)            = DefineSegments(Participant,Trial(i));
            Trial(i)            = ComputeKinematics(Trial(i),2,5); % Right lower limb kinematic chain
            Trial(i)            = ComputeKinematics(Trial(i),7,10); % Left lower limb kinematic chain
            Trial(i)            = ComputeKinematics(Trial(i),10,13); % Pelvis/lumbar/thorax/head
            Trial(i)            = ComputeKinematics(Trial(i),14,19); % Pelvis/lower lumbar/upper lumbar/lower thorax/upper thorax/head
            Trial(i).Joint(5)   = Trial(i).Joint(10); % Double pelvis/lumbar joint for indices coherence
            Trial(i).Segment(5) = Trial(i).Segment(10); % Double pelvis segment for indices coherence
        end
        
        % Process EMG signals
        Trial(i).EMG           = [];
        EMG                    = btkGetAnalogs(Trial(i).btk);
        Trial(i)               = InitialiseEMGSignals(Trial(i),EMG);
        fmethod.type           = 'butterBand4';
        fmethod.parameter      = [10 450];
        smethod.type           = 'butterLow2';
        smethod.parameter      = 3;
        [Calibration,Trial(i)] = ProcessEMGSignals(Calibration,Trial(i),0,fmethod,smethod,[]);
        clear EMG fmethod smethod;
        
        % Process forceplate signals
        Trial(i).GRF      = [];
        tGRF              = [];
        Trial(i).btk      = Correct_FP_C3D_Mokka(Trial(i).btk);
        tGRF              = btkGetForcePlatformWrenches(Trial(i).btk); % Required for C3D exportation only
        GRF               = btkGetGroundReactionWrenches(Trial(i).btk);
        GRFmeta           = btkGetMetaData(Trial(i).btk,'FORCE_PLATFORM');
        Trial(i)          = InitialiseGRFSignals(Trial(i),GRF,GRFmeta);
        fmethod.type      = 'threshold';
        fmethod.parameter = 35;
        smethod.type      = 'butterLow2';
        smethod.parameter = 50;
        [Trial(i),tGRF]   = ProcessGRFSignals(Session,Trial(i),GRF,tGRF,fmethod,smethod);
        clear GRF fmethod smethod;
                
        % Define additional events (for trials other than gait)
        % Crop raw files if needed to keep only wanted cycles
        if contains(Trial(i).type,'Perturbation_R_Shoulder')
            type      = 1;
            threshold = 135; % deg
            vec1      = Trial(i).Marker(58).Trajectory.smooth-Trial(i).Marker(54).Trajectory.smooth; % Vector RHAN-RSHO
            vec2      = Trial(i).Marker(8).Trajectory.smooth-Trial(i).Marker(54).Trajectory.smooth; % Vector RGTR-RSHO
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        
        elseif contains(Trial(i).type,'Perturbation_L_Shoulder')   
            type      = 1;
            threshold = 135; % deg
            vec1      = Trial(i).Marker(63).Trajectory.smooth-Trial(i).Marker(59).Trajectory.smooth; % Vector LHAN-lSHO
            vec2      = Trial(i).Marker(23).Trajectory.smooth-Trial(i).Marker(59).Trajectory.smooth; % Vector LGTR-LSHO
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        
        elseif contains(Trial(i).type,'Trunk_Forward')
            type      = 2;
            threshold = 45; % deg
            vec1      = Trial(i).Marker(41).Trajectory.smooth-...
                        (Trial(i).Marker(3).Trajectory.smooth+Trial(i).Marker(4).Trajectory.smooth)/2; % Vector C7-mean(RPSI,LPSI)
            vec2      = repmat([0 0 1],[size(vec1,1),1]); % Vector ICS_Z
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        
        elseif contains(Trial(i).type,'Trunk_Lateral')
            type      = 2;
            threshold = 22.5; % deg
            vec1      = Trial(i).Marker(41).Trajectory.smooth-...
                        (Trial(i).Marker(3).Trajectory.smooth+Trial(i).Marker(4).Trajectory.smooth)/2; % Vector C7-mean(RPSI,LPSI)
            vec2      = repmat([0 0 1],[size(vec1,1),1]); % Vector ICS_Z
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        
        elseif contains(Trial(i).type,'Trunk_Rotation')
            type      = 2;
            threshold = 22.5; % deg
            vec1      = Trial(i).Marker(54).Trajectory.smooth-Trial(i).Marker(59).Trajectory.smooth; % Vector RSHO-lSHO
            vec2      = Trial(i).Marker(15).Trajectory.smooth-Trial(i).Marker(30).Trajectory.smooth; % Vector RANK-LANK
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        
        elseif contains(Trial(i).type,'Weight_Constrained')
            type      = 2;
            threshold = 30; % deg
            vec1      = Trial(i).Marker(8).Trajectory.smooth-Trial(i).Marker(10).Trajectory.smooth; % Vector RGTR-RKNE
            vec2      = Trial(i).Marker(8).Trajectory.smooth-Trial(i).Marker(15).Trajectory.smooth; % Vector RANK-RKNE
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        
        elseif contains(Trial(i).type,'Weight_Unconstrained')
            type      = 2;
            threshold = 30; % deg
            vec1      = Trial(i).Marker(41).Trajectory.smooth-...
                        (Trial(i).Marker(3).Trajectory.smooth+Trial(i).Marker(4).Trajectory.smooth)/2; % Vector C7-mean(RPSI,LPSI)
            vec2      = repmat([0 0 1],[size(vec1,1),1]); % Vector ICS_Z
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        
        elseif contains(Trial(i).type,'S2S')
            type      = 3;
            threshold = 22.5; % deg
            vec1      = Trial(i).Marker(41).Trajectory.smooth-...
                        (Trial(i).Marker(3).Trajectory.smooth+Trial(i).Marker(4).Trajectory.smooth)/2; % Vector C7-mean(RPSI,LPSI)
            vec2      = repmat([0 0 1],[size(vec1,1),1]); % Vector ICS_Z
            Trial(i)  = DetectEvents(Trial(i),vec1,vec2,type,threshold);
            clear type threshold vec1 vec2;
        end
        
        % Cut data per cycle
        Trial(i) = CutCycles(Trial(i));

        % Store processed static data in a new C3D file
        ExportC3D(Trial(i),tGRF,Participant,Session,Folder);
        clear tGRF;
        
    end
end

% -------------------------------------------------------------------------
% COMPUTE BIOMARKERS
% -------------------------------------------------------------------------
disp('Compute biomarkers');

cd(Folder.biomarkers);
load('Biomarkers.mat');

% Store current participant and session in the list
Biomarker.participant = [Biomarker.participant [Participant.id,'_',Session.type(1:3)]];

% Set session indices
if contains(Participant.type,'Control')
    igroup = 1;
elseif contains(Participant.type,'Patient')
    igroup = 2;
end
iparticipant = str2num(Participant.id);
if contains(Session.type,'INI')
    isession = 1;
elseif contains(Session.type,'REL')
    isession = 2;
elseif contains(Session.type,'FWP')
    isession = 3;
end

% Compute and store participant/session biomarkers
% Biomarker dimensions : group x participant x session x side (dim 1 if
% central biomarker, dim 2 if right/left biomarker)

% BMo3
% d410 "Changing basic body position" Sit to stand	Pelvis/leg	Spatial/intensity	Hip sagittal angle (rom)
disp('  - BMo3');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        % Right side
        temp = [];
        for icycle = 1:size(Trial(i).Joint(4).Euler.rcycle,4)
            temp = [temp max(Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle),[],1) - ...
                         min(Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo3.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo3.units                                 = '°deg';
        clear temp;
        % Left side
        temp = [];
        for icycle = 1:size(Trial(i).Joint(9).Euler.rcycle,4)
            temp = [temp max(Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle),[],1) - ...
                         min(Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo3.value(igroup,iparticipant,isession,2) = rad2deg(mean(temp));
        Biomarker.BMo3.units                                 = '°deg';
        clear temp;
    end
end

% BMo4
% d410 "Changing basic body position" Stand to sit	Pelvis/leg	Spatial/intensity	Hip sagittal angle (rom)
disp('  - BMo4');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        % Right side
        temp = [];
        for icycle = 1:size(Trial(i).Joint(4).Euler.lcycle,4)
            temp = [temp max(Trial(i).Joint(4).Euler.lcycle(:,:,1,icycle),[],1) - ...
                         min(Trial(i).Joint(4).Euler.lcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo4.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo4.units                                 = '°deg';
        clear temp;
        % Left side
        temp = [];
        for icycle = 1:size(Trial(i).Joint(9).Euler.lcycle,4)
            temp = [temp max(Trial(i).Joint(9).Euler.lcycle(:,:,1,icycle),[],1) - ...
                         min(Trial(i).Joint(9).Euler.lcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo4.value(igroup,iparticipant,isession,2) = rad2deg(mean(temp));
        Biomarker.BMo4.units                                 = '°deg';
        clear temp;
    end
end

% BMo5
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar	Spatial/intensity	Lower lumbar sagittal angle (max)
disp('  - BMo5');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(15).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(15).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(15).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo5.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo5.units                                 = '°deg';
        clear temp;
    end
end

% BMo6
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar	Spatial/intensity	Lower lumbar sagittal angular velocity (max)
disp('  - BMo6');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(15).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(15).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(gradient(value),[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(15).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(gradient(value),[],1)];
                clear value;
            end
        end
        Biomarker.BMo6.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo6.units                                 = '°deg.frame-1';
        clear temp;
    end
end

% BMo9
% d410 "Changing basic body position" Trunk sagittal bending	Thorax	Spatial/intensity	Lower thorax sagittal angle (max)
disp('  - BMo9');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(17).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(17).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(17).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo9.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo9.units                                 = '°deg';
        clear temp;
    end
end

% BMo10
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar	Spatial/intensity	Lumbar contribution to thorax angle (rom)
disp('  - BMo10');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp1 = [];
        temp2 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(11).Euler.rcycle(:,1,2,icycle)-mean(pi/2-Trial(28).Segment(11).Euler.rcycle(1:10,1,2,icycle));
                temp1 = [temp1 max(value,[],1)];
                clear value;
                value = pi/2-Trial(i).Segment(12).Euler.rcycle(:,1,2,icycle)-mean(pi/2-Trial(28).Segment(12).Euler.rcycle(1:10,1,2,icycle));
                temp2 = [temp2 max(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(11).Euler.rcycle(:,1,1,icycle)-mean(pi/2-Trial(28).Segment(11).Euler.rcycle(1:10,1,1,icycle));
                temp1 = [temp1 max(value,[],1)];
                clear value;
                value = pi/2-Trial(i).Segment(12).Euler.rcycle(:,1,1,icycle)-mean(pi/2-Trial(28).Segment(12).Euler.rcycle(1:10,1,1,icycle));
                temp2 = [temp2 max(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo10.value(igroup,iparticipant,isession,1) = mean(temp1*100/temp2);
        Biomarker.BMo10.units                                 = '%';
        clear temp1 temp2;
    end
end

% BMo12
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar	Spatial/intensity	Lumbar sagittal angle (max)
disp('  - BMo12');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(11).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(11).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(11).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo12.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo12.units                                 = '°deg';
        clear temp;
    end
end

% BMo17
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar	Spatial/intensity	Lumbar sagittal angular velocity (max)
disp('  - BMo17');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(11).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(11).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(gradient(value),[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(11).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(gradient(value),[],1)];
                clear value;
            end
        end
        Biomarker.BMo17.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo17.units                                 = '°deg.frame-1';
        clear temp;
    end
end

% BMo23
% d410 "Changing basic body position" Sit to stand	Lumbar/leg	Spatial/intensity	Lumbar/hip ratio of sagittal angle (rom)
disp('  - BMo23');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        temp1 = [];
        temp2 = [];
        temp3 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            temp1 = [temp1 max(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle),[],1) - ...
                           min(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle),[],1)];
            temp2 = [temp2 max(Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle),[],1) - ...
                           min(Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle),[],1)];
            temp3 = [temp3 max(Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle),[],1) - ...
                           min(Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo23.value(igroup,iparticipant,isession,1) = mean(temp1/temp2);
        Biomarker.BMo23.value(igroup,iparticipant,isession,2) = mean(temp1/temp3);
        Biomarker.BMo23.units                                 = 'ratio';
        clear temp1 temp2;
    end
end

% BMo24
% d410 "Changing basic body position" Stand to sit	Lumbar/leg	Spatial/intensity	Lumbar/hip ratio of sagittal angle (rom)
disp('  - BMo24');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        temp1 = [];
        temp2 = [];
        temp3 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.lcycle,4)
            temp1 = [temp1 max(Trial(i).Joint(10).Euler.lcycle(:,:,1,icycle),[],1) - ...
                           min(Trial(i).Joint(10).Euler.lcycle(:,:,1,icycle),[],1)];
            temp2 = [temp2 max(Trial(i).Joint(4).Euler.lcycle(:,:,1,icycle),[],1) - ...
                           min(Trial(i).Joint(4).Euler.lcycle(:,:,1,icycle),[],1)];
            temp3 = [temp3 max(Trial(i).Joint(9).Euler.lcycle(:,:,1,icycle),[],1) - ...
                           min(Trial(i).Joint(9).Euler.lcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo24.value(igroup,iparticipant,isession,1) = mean(temp1/temp2);
        Biomarker.BMo24.value(igroup,iparticipant,isession,2) = mean(temp1/temp3);
        Biomarker.BMo24.units                                 = 'ratio';
        clear temp1 temp2;
    end
end

% BMo25
% d410 "Changing basic body position" Sit to stand	Lumbar/leg	Coordination	Lumbar/hip relative phase difference (max)
disp('  - BMo25');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        phi1 = [];
        phi2 = [];
        phi3 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            phi1 = [phi1 atan(gradient(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle))];
            phi2 = [phi2 atan(gradient(Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle))];
            phi3 = [phi3 atan(gradient(Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle))];
        end
        relativePhase1 = mean(phi2-phi1,2);
        relativePhase2 = mean(phi3-phi1,2);
        Biomarker.BMo25.value(igroup,iparticipant,isession,1) = rad2deg(max(relativePhase1));
        Biomarker.BMo25.value(igroup,iparticipant,isession,2) = rad2deg(max(relativePhase2));
        Biomarker.BMo25.units                                 = '°deg';
        clear phi1 phi2 phi3 relativePhase1 relativePhase2;        
    end
end

% BMo26
% d410 "Changing basic body position" Sit to stand	Lumbar/leg	Coordination	Lumbar/hip relative phase difference (mean)
disp('  - BMo26');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        phi1 = [];
        phi2 = [];
        phi3 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            phi1 = [phi1 atan(gradient(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle))];
            phi2 = [phi2 atan(gradient(Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle))];
            phi3 = [phi3 atan(gradient(Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle))];
        end
        relativePhase1 = mean(phi2-phi1,2);
        relativePhase2 = mean(phi3-phi1,2);
        Biomarker.BMo26.value(igroup,iparticipant,isession,1) = rad2deg(mean(relativePhase1));
        Biomarker.BMo26.value(igroup,iparticipant,isession,2) = rad2deg(mean(relativePhase2));
        Biomarker.BMo26.units                                 = '°deg';
        clear phi1 phi2 phi3 relativePhase1 relativePhase2;        
    end
end

% BMo27
% d410 "Changing basic body position" Stand to sit	Lumbar/leg	Coordination	Lumbar/hip relative phase difference (mean)
disp('  - BMo27');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        phi1 = [];
        phi2 = [];
        phi3 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.lcycle,4)
            phi1 = [phi1 atan(gradient(Trial(i).Joint(10).Euler.lcycle(:,:,1,icycle))./Trial(i).Joint(10).Euler.lcycle(:,:,1,icycle))];
            phi2 = [phi2 atan(gradient(Trial(i).Joint(4).Euler.lcycle(:,:,1,icycle))./Trial(i).Joint(4).Euler.lcycle(:,:,1,icycle))];
            phi3 = [phi3 atan(gradient(Trial(i).Joint(9).Euler.lcycle(:,:,1,icycle))./Trial(i).Joint(9).Euler.lcycle(:,:,1,icycle))];
        end
        relativePhase1 = mean(phi2-phi1,2);
        relativePhase2 = mean(phi3-phi1,2);
        Biomarker.BMo27.value(igroup,iparticipant,isession,1) = rad2deg(mean(relativePhase1));
        Biomarker.BMo27.value(igroup,iparticipant,isession,2) = rad2deg(mean(relativePhase2));
        Biomarker.BMo27.units                                 = '°deg';
        clear phi1 phi2 phi3 relativePhase1 relativePhase2;        
    end
end

% BMo28
% d410 "Changing basic body position" Sit to stand	Lumbar/leg	Coordination	Lumbar/hip relative phase difference (min)
disp('  - BMo28');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        phi1 = [];
        phi2 = [];
        phi3 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            phi1 = [phi1 atan(gradient(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle))];
            phi2 = [phi2 atan(gradient(Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(4).Euler.rcycle(:,:,1,icycle))];
            phi3 = [phi3 atan(gradient(Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle))./Trial(i).Joint(9).Euler.rcycle(:,:,1,icycle))];
        end
        relativePhase1 = mean(phi2-phi1,2);
        relativePhase2 = mean(phi3-phi1,2);
        Biomarker.BMo28.value(igroup,iparticipant,isession,1) = rad2deg(min(relativePhase1));
        Biomarker.BMo28.value(igroup,iparticipant,isession,2) = rad2deg(min(relativePhase2));
        Biomarker.BMo28.units                                 = '°deg';
        clear phi1 phi2 phi3 relativePhase1 relativePhase2;        
    end
end

% BMo29
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar/pelvis	Coordination	Lumbar/pelvis absolute relative phase (mean)
disp('  - BMo29');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        phi1 = [];
        phi2 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                ang1  = pi/2-Trial(i).Segment(10).Euler.rcycle(:,:,2,icycle);
                nang1 = 2*(ang1-min(ang1))-(max(ang1)-min(ang1))-1;
                vel1  = gradient(ang1);
                nvel1 = vel1/max(vel1);
                phi1  = [phi1 atan(nvel1./nang1)];
                ang2  = pi/2-Trial(i).Segment(11).Euler.rcycle(:,:,2,icycle);
                nang2 = 2*(ang2-min(ang2))-(max(ang2)-min(ang2))-1;
                vel2  = gradient(ang2);
                nvel2 = vel2/max(vel2);
                phi2  = [phi2 atan(nvel2./nang2)];
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                ang1  = pi/2-Trial(i).Segment(10).Euler.rcycle(:,:,1,icycle);
                nang1 = 2*(ang1-min(ang1))-(max(ang1)-min(ang1))-1;
                vel1  = gradient(ang1);
                nvel1 = vel1/max(vel1);
                phi1  = [phi1 atan(nvel1./nang1)];
                ang2  = pi/2-Trial(i).Segment(11).Euler.rcycle(:,:,1,icycle);
                nang2 = 2*(ang2-min(ang2))-(max(ang2)-min(ang2))-1;
                vel2  = gradient(ang2);
                nvel2 = vel2/max(vel2);
                phi2  = [phi2 atan(nvel2./nang2)];
            end
        end
        relativePhase = mean(phi1-phi2,2);
        MARP          = sum(abs(relativePhase),1)/101;
        Biomarker.BMo29.value(igroup,iparticipant,isession,1) = rad2deg(MARP);
        Biomarker.BMo29.units                                 = '°deg';
        clear ang1 nang1 vel1 nvel1 ang2 nang2 vel2 nvel2 phi1 phi2 relativePhase MARP;  
    end
end

% BMo30
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar/pelvis	Coordination	Lumbar/pelvis deviation phase (mean)
disp('  - BMo30');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        phi1 = [];
        phi2 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                ang1  = pi/2-Trial(i).Segment(10).Euler.rcycle(:,:,2,icycle);
                nang1 = 2*(ang1-min(ang1))-(max(ang1)-min(ang1))-1;
                vel1  = gradient(ang1);
                nvel1 = vel1/max(vel1);
                phi1  = [phi1 atan(nvel1./nang1)];
                ang2  = pi/2-Trial(i).Segment(11).Euler.rcycle(:,:,2,icycle);
                nang2 = 2*(ang2-min(ang2))-(max(ang2)-min(ang2))-1;
                vel2  = gradient(ang2);
                nvel2 = vel2/max(vel2);
                phi2  = [phi2 atan(nvel2./nang2)];
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                ang1  = pi/2-Trial(i).Segment(10).Euler.rcycle(:,:,1,icycle);
                nang1 = 2*(ang1-min(ang1))-(max(ang1)-min(ang1))-1;
                vel1  = gradient(ang1);
                nvel1 = vel1/max(vel1);
                phi1  = [phi1 atan(nvel1./nang1)];
                ang2  = pi/2-Trial(i).Segment(11).Euler.rcycle(:,:,1,icycle);
                nang2 = 2*(ang2-min(ang2))-(max(ang2)-min(ang2))-1;
                vel2  = gradient(ang2);
                nvel2 = vel2/max(vel2);
                phi2  = [phi2 atan(nvel2./nang2)];
            end
        end
        relativePhase = phi1-phi2;
        DP            = sum(std(relativePhase,0,2),1)/101;
        Biomarker.BMo30.value(igroup,iparticipant,isession,1) = rad2deg(DP);
        Biomarker.BMo30.units                                 = '°deg';
        clear ang1 nang1 vel1 nvel1 ang2 nang2 vel2 nvel2 phi1 phi2 relativePhase DP;  
    end
end

% BMo33
% d410 "Changing basic body position" Sit to stand	Lumbar/pelvis	Spatial/intensity	Lumbopelvic sagittal angle (rom)
disp('  - BMo33');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        temp = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            temp = [temp max(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle),[],1) - ...
                         min(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo33.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo33.units                                 = '°deg';
        clear temp;
    end
end

% BMo34
% d410 "Changing basic body position" Stand to sit	Lumbar/pelvis	Spatial/intensity	Lumbopelvic sagittal angle (rom)
disp('  - BMo34');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'S2S_Unconstrained')
        temp = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.lcycle,4)
            temp = [temp max(Trial(i).Joint(10).Euler.lcycle(:,:,1,icycle),[],1) - ...
                         min(Trial(i).Joint(10).Euler.lcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo34.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo34.units                                 = '°deg';
        clear temp;
    end
end

% BMo37
% d410 "Changing basic body position" Trunk sagittal bending	Pelvis	Spatial/intensity	Pelvis sagittal angle (max)
disp('  - BMo37');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(10).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(value,[],1)-min(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(10).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(value,[],1)-min(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo37.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo37.units                                 = '°deg';
        clear temp;
    end
end

% BMo42
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar/pelvis	Coordination	Lumbar/pelvis deviation phase (mean)
disp('  - BMo42');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        phi1 = [];
        phi2 = [];
        phi3 = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                ang1  = pi/2-Trial(i).Segment(4).Euler.rcycle(:,:,2,icycle);
                nang1 = 2*(ang1-min(ang1))-(max(ang1)-min(ang1))-1;
                vel1  = gradient(ang1);
                nvel1 = vel1/max(vel1);
                phi1  = [phi1 atan(nvel1./nang1)];
                ang2  = pi/2-Trial(i).Segment(9).Euler.rcycle(:,:,2,icycle);
                nang2 = 2*(ang2-min(ang2))-(max(ang2)-min(ang2))-1;
                vel2  = gradient(ang2);
                nvel2 = vel1/max(vel2);
                phi2  = [phi2 atan(nvel2./nang2)];
                ang3  = pi/2-Trial(i).Segment(10).Euler.rcycle(:,:,2,icycle);
                nang3 = 2*(ang2-min(ang3))-(max(ang3)-min(ang3))-1;
                vel3  = gradient(ang3);
                nvel3 = vel2/max(vel3);
                phi3  = [phi3 atan(nvel3./nang3)];
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                ang1  = pi/2-Trial(i).Segment(4).Euler.rcycle(:,:,1,icycle);
                nang1 = 2*(ang1-min(ang1))-(max(ang1)-min(ang1))-1;
                vel1  = gradient(ang1);
                nvel1 = vel1/max(vel1);
                phi1  = [phi1 atan(nvel1./nang1)];
                ang2  = pi/2-Trial(i).Segment(9).Euler.rcycle(:,:,1,icycle);
                nang2 = 2*(ang2-min(ang2))-(max(ang2)-min(ang2))-1;
                vel2  = gradient(ang2);
                nvel2 = vel1/max(vel2);
                phi2  = [phi2 atan(nvel2./nang2)];
                ang3  = pi/2-Trial(i).Segment(10).Euler.rcycle(:,:,1,icycle);
                nang3 = 2*(ang2-min(ang3))-(max(ang3)-min(ang3))-1;
                vel3  = gradient(ang3);
                nvel3 = vel2/max(vel3);
                phi3  = [phi3 atan(nvel3./nang3)];
            end
        end
        relativePhase1 = phi1-phi3;
        DP1            = sum(std(relativePhase1,0,2),1)/101;
        relativePhase2 = phi2-phi3;
        DP2            = sum(std(relativePhase2,0,2),1)/101;
        Biomarker.BMo42.value(igroup,iparticipant,isession,1) = rad2deg(DP1);
        Biomarker.BMo42.value(igroup,iparticipant,isession,2) = rad2deg(DP2);
        Biomarker.BMo42.units                                 = '°deg';
        clear ang1 nang1 vel1 nvel1 phi1 ang2 nang2 vel2 nvel2 phi2 ang3 nang3 vel3 nvel3 phi3 relativePhase1 relativePhase2 DP1 DP2;  
    end
end

% BMo43
% d410 "Changing basic body position" Trunk rotation	Thorax	Spatial/intensity	Scapular belt transversal angle (max)
disp('  - BMo43');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Rotation')
        temp = [];
        for icycle = 1:size(Trial(i).Marker(54).Trajectory.rcycle,3)
            vec1 = [];
            vec2 = [];
            ang  = [];
            vec1 = Trial(i).Marker(54).Trajectory.rcycle(:,1:2,icycle) - Trial(i).Marker(59).Trajectory.rcycle(:,1:2,icycle);
            vec1 = [vec1 zeros(101,1)];
            vec2 = Trial(i).Marker(1).Trajectory.rcycle(:,1:2,icycle) - Trial(i).Marker(6).Trajectory.rcycle(:,1:2,icycle);
            vec2 = [vec2 zeros(101,1)];
            for t = 1:101
                ang(t)  = atan2(norm(cross(vec1(t,:),vec2(t,:))),dot(vec1(t,:),vec2(t,:)));
            end
            temp = [temp max(ang)];
        end
        Biomarker.BMo43.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo43.units                                 = '°deg';
        clear temp vec1 vec2 ang;
    end
end

% BMo44
% d410 "Changing basic body position" Trunk sagittal bending	Thorax/pelvis	Spatial/intensity	Thoracopelvic sagittal angle (max)
disp('  - BMo44');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.lcycle,4)
            temp = [temp max(abs(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle)) + ...
                             abs(Trial(i).Joint(11).Euler.rcycle(:,:,1,icycle)))];
        end
        Biomarker.BMo44.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo44.units                                 = '°deg';
        clear temp;
    end
end

% BMo49
% d410 "Changing basic body position" Trunk sagittal bending	Thorax	Spatial/intensity	Thorax sagittal angle (rom)
disp('  - BMo49');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(12).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(12).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(value,[],1)-min(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(12).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(value,[],1)-min(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo49.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo49.units                                 = '°deg';
        clear temp;
    end
end

% BMo57
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar	Spatial/intensity	Upper lumbar sagittal angle (max)
disp('  - BMo57');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(16).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(16).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(16).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo57.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo57.units                                 = '°deg';
        clear temp;
    end
end

% BMo58
% d410 "Changing basic body position" Trunk sagittal bending	Lumbar	Spatial/intensity	Upper lumbar sagittal angular velocity (max)
disp('  - BMo58');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(16).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(16).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(gradient(value),[],1)-min(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(16).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(gradient(value),[],1)-min(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo58.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo58.units                                 = '°deg.frame-1';
        clear temp;
    end
end

% BMo59
% d410 "Changing basic body position" Trunk sagittal bending	Thorax	Spatial/intensity	Upper thorax sagittal angle (max)
disp('  - BMo59');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Trunk_Forward')
        temp = [];
        for icycle = 1:size(Trial(i).Segment(18).Euler.rcycle,4)
            if abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,1),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,1),1))>0.1 % Y forward
                value = pi/2-Trial(i).Segment(18).Euler.rcycle(:,1,2,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            elseif abs(mean(Trial(i).Marker(1).Trajectory.smooth(:,2),1)-mean(Trial(i).Marker(6).Trajectory.smooth(:,2),1))>0.1 % X forward
                value = pi/2-Trial(i).Segment(18).Euler.rcycle(:,1,1,icycle);
                temp  = [temp max(value,[],1)];
                clear value;
            end
        end
        Biomarker.BMo59.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo59.units                                 = '°deg';
        clear temp;
    end
end

% BMo80
% d430 "Lifting and carrying objects" Trunk sagittal bending	Lumbar/pelvis	Spatial/intensity	Lumbopelvic sagittal angle (max)
disp('  - BMo80');
for i = 1:size(Trial,2)
    if contains(Trial(i).type,'Weight_Unconstrained')
        temp = [];
        for icycle = 1:size(Trial(i).Joint(10).Euler.rcycle,4)
            temp = [temp max(Trial(i).Joint(10).Euler.rcycle(:,:,1,icycle),[],1)];
        end
        Biomarker.BMo80.value(igroup,iparticipant,isession,1) = rad2deg(mean(temp));
        Biomarker.BMo80.units                                 = '°deg';
        clear temp;
    end
end