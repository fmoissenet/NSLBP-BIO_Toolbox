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
              'sMVC','Test'};

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
    
    % Events
    Static(i).Event = [];
    
    % Marker trajectories
    Marker           = btkGetMarkers(Static(i).btk);
    Static(i).Marker = [];
    Static(i)        = InitialiseMarkerTrajectories(Static(i),Marker);
    Static(i)        = ProcessMarkerTrajectories([],Static(i));
    clear Marker;
    
    % EMG signals
    Static(i).EMG = [];
    
    % Store processed static data in a new C3D file
    ExportC3D(Static(i),Participant,Session,Folder);
end

% EMG calibration data
disp('Pre-process EMG calibration data');
Calibration(1).type = 'EMG_calibration';
Calibration(1).EMG  = [];
for i = 1:size(Trial,2)
    
    if ~isempty(strfind(Trial(i).type,'Endurance')) || ...
       ~isempty(strfind(Trial(i).type,'sMVC')) 
        
        disp(['  - ',Trial(i).file]);

        % Events
        Event          = btkGetEvents(Trial(i).btk);
        Trial(i).Event = [];
        Trial(i)       = InitialiseEvents(Trial(i),Event);
        clear Event;   
    
        % Marker trajectories
        Trial(i).Marker = [];  
        
        % EMG signals
        EMG                    = btkGetAnalogs(Trial(i).btk);
        Trial(i).EMG           = [];
        Trial(i)               = InitialiseEMGSignals(Trial(i),EMG);
        fmethod.type           = 'butterBand4';
        fmethod.parameter      = [10 450];
        smethod.type           = 'butterLow2';
        smethod.parameter      = 3;
        nmethod.type           = 'sMVC';
        [Calibration,Trial(i)] = ProcessEMGSignals(Calibration,Trial(i),1,fmethod,smethod,nmethod);
        clear EMG;
        
        % Store processed static data in a new C3D file
        ExportC3D(Trial(i),Participant,Session,Folder);        
   
    end
end

% Trial data
disp('Pre-process trial data');
for i = [12,28]%1:size(Trial,2)
    
    if isempty(strfind(Trial(i).type,'sMVC')) % Endurance tasks considered as Trial here
        
        disp(['  - ',Trial(i).file]);

        % Events
        Event          = btkGetEvents(Trial(i).btk);
        Trial(i).Event = [];
        Trial(i)       = InitialiseEvents(Trial(i),Event);
        clear Event;   

        % Marker trajectories   
        Marker               = btkGetMarkers(Trial(i).btk);
        Trial(i).Marker      = [];
        Trial(i)             = InitialiseMarkerTrajectories(Trial(i),Marker);        
        fmethod.type         = 'intercor';
        fmethod.gapThreshold = [];
        smethod.type         = 'movmean';
        smethod.parameter    = 10;        
        Trial(i)             = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);        
        clear Marker;
                
        % EMG signals
        EMG                    = btkGetAnalogs(Trial(i).btk);
        Trial(i).EMG           = [];
        Trial(i)               = InitialiseEMGSignals(Trial(i),EMG);
        fmethod.type           = 'butterBand4';
        fmethod.parameter      = [10 450];
        smethod.type           = 'butterLow2';
        smethod.parameter      = 3;
        [Calibration,Trial(i)] = ProcessEMGSignals(Calibration,Trial(i),0,fmethod,smethod,[]);
        clear EMG;

        % Store processed static data in a new C3D file
        ExportC3D(Trial(i),Participant,Session,Folder);
        
    end
end