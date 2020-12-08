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
% Description  : This routine aims to detect foot / force plate contacts.
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

function detectFootForceplate(Marker,Grf,Event)

if ~isempty(Event) % Events only used for gait trials
    
    % Right foot
    for i = 1:size(Event.RHS,2)
        % Get foot markers position at foot strike and foot off
        RFootMarkersStart = [Marker.RHEE(Event.RHS(i),:); ...
                             Marker.RFMH(Event.RHS(i),:); ...
                             Marker.RVMH(Event.RHS(i),:)];
        if Event.RTO(i) > Event.RHS(i)
            RFootMarkersStop  = [Marker.RHEE(Event.RTO(i),:); ...
                                 Marker.RFMH(Event.RTO(i),:); ...
                                 Marker.RVMH(Event.RTO(i),:)];
        else
            RFootMarkersStop  = [Marker.RHEE(Event.RTO(i+1),:); ...
                                 Marker.RFMH(Event.RTO(i+1),:); ...
                                 Marker.RVMH(Event.RTO(i+1),:)];
        end
        
        % Forceplace 1: Check if all foot markers (+threshold) included in 
        % the forceplate corners
        if RFootMarkersStart(1,1) > Grf(1).corner(1,1) && RFootMarkersStart(1,1) > Grf(1).corner(2,1)
    end
    
end