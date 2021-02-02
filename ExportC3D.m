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
% Description  : This routine aims to export C3D files with updated data.
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

function ExportC3D(Trial,tGRF,Participant,Session,Folder)

% Set new C3D file
btkFile = btkNewAcquisition();
btkSetFrequency(btkFile,Trial.fmarker);
btkSetFrameNumber(btkFile,Trial.n1);
btkSetPointsUnit(btkFile,'marker','m');
btkSetAnalogSampleNumberPerFrame(btkFile,10);

% Append events
if ~isempty(Trial.Event)
    for i = 1:size(Trial.Event,2)
        for j = 1:size(Trial.Event(i).value,2)
            Event = Trial.Event(i).value(1,j)/Trial.fmarker;
            btkAppendEvent(btkFile,Trial.Event(i).label,Event,'');
            clear Event;
        end
    end
end

% Append marker trajectories
if ~isempty(Trial.Marker)
    for i = 1:size(Trial.Marker,2)
        if ~isempty(Trial.Marker(i).Trajectory.smooth)
            btkAppendPoint(btkFile,'marker',Trial.Marker(i).label,Trial.Marker(i).Trajectory.smooth);
        else
            btkAppendPoint(btkFile,'marker',Trial.Marker(i).label,zeros(Trial.n1,3));
        end
    end
end

% Append EMG signals
if ~isempty(Trial.EMG)
    for i = 1:size(Trial.EMG,2)
        if ~isempty(Trial.EMG(i).Signal.norm)
            btkAppendAnalog(btkFile,Trial.EMG(i).label,Trial.EMG(i).Signal.norm,'EMG signal (normalised)');
        else
            btkAppendAnalog(btkFile,Trial.EMG(i).label,Trial.EMG(i).Signal.smooth,'EMG signal (mV)');
        end
    end
end

% Append GRF signals
if ~isempty(Trial.GRF)
    GRF     = btkGetGroundReactionWrenches(Trial.btk);
    GRFmeta = btkGetMetaData(Trial.btk,'FORCE_PLATFORM');
    for i = 1:size(GRF,1)
        GRF(i).corners = GRFmeta.children.CORNERS.info.values(:,:,i)*1e-3;
        GRF(i).origin  = GRFmeta.children.ORIGIN.info.values(:,i)*1e-3;
        Trial.GRF(i).Signal.M.smooth = Trial.GRF(i).Signal.M.smooth;
        if ~isempty(Trial.GRF(i).Signal.F.smooth)
            btkAppendForcePlatformType2(btkFile,tGRF(i).F,...
                                        tGRF(i).M,...
                                        GRF(i).corners',GRF(i).origin',[0,0,0]);
        else
            btkAppendForcePlatformType2(btkFile,zeros(size(GRF(1).F)),...
                                        zeros(size(GRF(1).M)),...
                                        GRF(i).corners',GRF(i).origin',[0,0,0]);
        end
    end
end

% Append participant metadata
nData = 9;
info.format = 'Integer';
info.values = nData;
btkAppendMetaData(btkFile,'PARTICIPANT','USED',info);
clear info;
info.format = 'Char';
info.dimensions = ['1x',nData];
info.values(1:nData) = {'gender' 'inclusionAge' 'pelvisWidth' ...
                        'RLegLength' 'LLegLength' ...
                        'RKneeWidth' 'LKneeWidth' ...
                        'RAnkleWidth' 'LAnkleWidth'};
btkAppendMetaData(btkFile,'PARTICIPANT','LABELS',info);
clear info;
info.format = 'Char';
info.dimensions = ['1x',nData];
info.values(1:nData) = {'adimensioned (0: female, 1: male)' 'years' 'm' ...
                        'm' 'm' ...
                        'm' 'm' ...
                        'm' 'm'};
btkAppendMetaData(btkFile,'PARTICIPANT','UNITS',info);
clear info;
info.format     = 'Real';
info.dimensions = ['1x',nData];
if strcmp(Participant.gender,'Male')
    info.values(2) = 1;
elseif strcmp(Participant.gender,'Female')
    info.values(2) = 0;
end
info.values(3)  = Participant.inclusionAge;
info.values(4)  = Participant.pelvisWidth*1e-2;
info.values(5)  = Participant.RLegLength*1e-2;
info.values(6)  = Participant.LLegLength*1e-2;
info.values(7)  = Participant.RKneeWidth*1e-2;
info.values(8)  = Participant.LKneeWidth*1e-2;
info.values(9)  = Participant.RAnkleWidth*1e-2;
info.values(10) = Participant.LAnkleWidth*1e-2;
btkAppendMetaData(btkFile,'PARTICIPANT','VALUES',info);

% Append session metadata
nData                = 5;
info.format          = 'Integer';
info.values          = nData;
btkAppendMetaData(btkFile,'SESSION','USED',info);
clear info;
info.format          = 'Char';
info.dimensions      = ['1x',nData];
info.values(1:nData) = {'date' 'type' 'examiner' ...
                        'participantHeight' 'participantWeight'};
btkAppendMetaData(btkFile,'SESSION','LABELS',info);
clear info;
info.format          = 'Char';
info.dimensions      = ['1x',nData];
info.values(1:nData) = {'DD-MM-YYYY' 'XXX_session' 'initials' 'm' 'kg'};
btkAppendMetaData(btkFile,'SESSION','UNITS',info);
clear info;
info.format          = 'Char';
info.dimensions      = ['1x',nData];
info.values(1:nData) = {Session.date Session.type Session.examiner ...
                        num2str(Session.participantHeight*1e-2) num2str(Session.participantWeight)};
btkAppendMetaData(btkFile,'SESSION','VALUES',info);

% Export C3D file
cd(Folder.export);
btkWriteAcquisition(btkFile,[regexprep(Trial.file,'.c3d',''),'_processed.c3d']);