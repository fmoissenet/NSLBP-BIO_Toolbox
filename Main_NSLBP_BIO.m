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
% Inputs       : To be defined
% Outputs      : To be defined
% -------------------------------------------------------------------------
% Dependencies : - Biomechanical Toolkit (BTK): https://github.com/Biomechanical-ToolKit/BTKCore
%                - Soder: I. Soederqvist and P.A. Wedin (1993) Determining the movement of the skeleton
%                         using well-configured markers. J. Biomech. 26:1473-1477
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
addpath(genpath(Folder.dependencies)); % Adds the specified folders to the top of the search path for the current session

% -------------------------------------------------------------------------
% SET PARTICIPANT
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

% List all the trial types
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
k1 = 1;
k2 = 1;
for i = 1:size(c3dFiles)
    disp(['  - ',c3dFiles(i).name]);
    for j = 1:size(trialTypes,2)
        if strfind(c3dFiles(i).name,trialTypes{j}) == 1
            if strcmp(trialTypes{j},'Static')
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
% PRE-PROCESS STATIC DATA
% -------------------------------------------------------------------------
disp('Pre-process static data');
for i = 1:size(Static,2)
    disp(['  - ',Static(i).file]);
    
    % Marker trajectories
    Marker               = btkGetMarkers(Static(i).btk);
    Static(i).Marker     = [];
    Static(i)            = InitialiseMarkerset(Static(i),Marker);
    Static(i)            = ProcessMarkerTrajectories([],Static(i));
    clear Marker;
    
    % Store processed static data in a new C3D file
%     ExportC3D(Static(i),Participant,Session,Folder);
end

% -------------------------------------------------------------------------
% PRE-PROCESS TRIAL DATA
% -------------------------------------------------------------------------
disp('Pre-process trial data');
for i = 1%1:size(Trial,2)
    disp(['  - ',Trial(i).file]);
    
    % Marker trajectories
    if isempty(strfind(Trial(i).type,'Endurance')) && ...
       isempty(strfind(Trial(i).type,'sMVC'))   
        Marker               = btkGetMarkers(Trial(i).btk);
        Trial(i).Marker      = [];
        Trial(i)             = InitialiseMarkerset(Trial(i),Marker);
%         fmethod.type         = 'none';
%         fmethod.gapThreshold = [];
%         smethod.type         = 'none';
%         smethod.parameter    = [];        
%         Trial(i)             = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);
        
        fmethod.type         = 'none';
        fmethod.gapThreshold = [];
        smethod.type         = 'none';
        smethod.parameter    = [];        
        TrialREF(i)          = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);

        % Create gaps on RASI
        gapSize  = 20; % frames
        gapSpace = 200; % frames
        for t = 1:gapSpace:Trial(i).n1-10
            Trial(i).Marker(1).Trajectory.raw(t:t+gapSize-1,:) = zeros(size(Trial(i).Marker(1).Trajectory.raw(t:t+gapSize-1,:)));
        end
        
        fmethod.type         = 'none';
        fmethod.gapThreshold = [];
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial00(i)          = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);
        
        fmethod.type         = 'linear';
        fmethod.gapThreshold = 50;
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial10(i)           = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);

        fmethod.type         = 'spline';
        fmethod.gapThreshold = 50;
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial20(i)           = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);

        fmethod.type         = 'pchip';
        fmethod.gapThreshold = 50;
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial30(i)           = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);

        fmethod.type         = 'makima';
        fmethod.gapThreshold = 50;
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial40(i)           = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);

        fmethod.type         = 'intercor';
        fmethod.gapThreshold = [];
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial50(i)           = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);

        fmethod.type         = 'rigid';
        fmethod.gapThreshold = [];
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial60(i)           = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);

        fmethod.type         = 'kalman';
        fmethod.gapThreshold = [];
        smethod.type         = 'none';
        smethod.parameter    = [];        
        Trial70(i)           = ProcessMarkerTrajectories(Static,Trial(i),fmethod,smethod);
        
        figure; hold on;
        subplot(1,3,1); hold on; title('X');
        plot(TrialREF(i).Marker(1).Trajectory.smooth(:,1),'black');
        plot(Trial00(i).Marker(1).Trajectory.smooth(:,1),'red');
        plot(Trial10(i).Marker(1).Trajectory.smooth(:,1),'green');
        plot(Trial20(i).Marker(1).Trajectory.smooth(:,1),'blue');
        plot(Trial30(i).Marker(1).Trajectory.smooth(:,1),'yellow');
        plot(Trial40(i).Marker(1).Trajectory.smooth(:,1),'cyan');
        plot(Trial50(i).Marker(1).Trajectory.smooth(:,1),'magenta');
        plot(Trial60(i).Marker(1).Trajectory.smooth(:,1),'--r');
        plot(Trial70(i).Marker(1).Trajectory.smooth(:,1),'--g');
        subplot(1,3,2); hold on; title('Y');
        plot(TrialREF(i).Marker(1).Trajectory.smooth(:,2),'black');
        plot(Trial00(i).Marker(1).Trajectory.smooth(:,2),'red');
        plot(Trial10(i).Marker(1).Trajectory.smooth(:,2),'green');
        plot(Trial20(i).Marker(1).Trajectory.smooth(:,2),'blue');
        plot(Trial30(i).Marker(1).Trajectory.smooth(:,2),'yellow');
        plot(Trial40(i).Marker(1).Trajectory.smooth(:,2),'cyan');
        plot(Trial50(i).Marker(1).Trajectory.smooth(:,2),'magenta');
        plot(Trial60(i).Marker(1).Trajectory.smooth(:,2),'--r');
        plot(Trial70(i).Marker(1).Trajectory.smooth(:,2),'--g');
        subplot(1,3,3); hold on; title('Z');
        plot(TrialREF(i).Marker(1).Trajectory.smooth(:,3),'black');
        plot(Trial00(i).Marker(1).Trajectory.smooth(:,3),'red');
        plot(Trial10(i).Marker(1).Trajectory.smooth(:,3),'green');
        plot(Trial20(i).Marker(1).Trajectory.smooth(:,3),'blue');
        plot(Trial30(i).Marker(1).Trajectory.smooth(:,3),'yellow');
        plot(Trial40(i).Marker(1).Trajectory.smooth(:,3),'cyan');
        plot(Trial50(i).Marker(1).Trajectory.smooth(:,3),'magenta');
        plot(Trial60(i).Marker(1).Trajectory.smooth(:,3),'--r');
        plot(Trial70(i).Marker(1).Trajectory.smooth(:,3),'--g');

        clear Marker;
    end

    % EMG signals
    if ~isempty(strfind(Trial(i).type,'Endurance')) || ...
       ~isempty(strfind(Trial(i).type,'sMVC')) 
        Trial(i).Marker = [];
    end
        
    % Store processed static data in a new C3D file
%     ExportC3D(Trial(i),Participant,Session,Folder);

end