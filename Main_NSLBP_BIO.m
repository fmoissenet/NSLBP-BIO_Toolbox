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
addpath(Folder.toolbox);
addpath(genpath(Folder.dependencies));

% -------------------------------------------------------------------------
% DEFINE PARTICIPANT
% -------------------------------------------------------------------------
disp('Set participant parameters');
Participant.gender       = 'Male'; % (0: woman, 1: man)
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