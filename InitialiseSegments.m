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
% Dependencies : None
% -------------------------------------------------------------------------
% This work is licensed under the Creative Commons Attribution - 
% NonCommercial 4.0 International License. To view a copy of this license, 
% visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% -------------------------------------------------------------------------

function Trial = InitialiseSegments(Trial)

segmentLabels = {'Right forceplate','Right foot','Right tibia','Right femur','Pelvis',...
                 'Left forceplate','Left foot','Left tibia','Left femur','Pelvis',...
                 'Lumbar','Thorax','Head',...
                 'Pelvis',...
                 'Lower lumbar','Upper lumbar','Lower thorax','Upper thorax','Head'};
             
for i = 1:length(segmentLabels)
    Trial.Segment(i).label        = segmentLabels{i};
    Trial.Segment(i).rM.smooth    = [];
    Trial.Segment(i).rM.rcycle    = [];
    Trial.Segment(i).rM.lcycle    = [];
    Trial.Segment(i).rM.units     = 'm';
    Trial.Segment(i).Q.smooth     = [];
    Trial.Segment(i).Q.rcycle     = [];
    Trial.Segment(i).Q.lcycle     = [];
    Trial.Segment(i).T.smooth     = [];
    Trial.Segment(i).T.rcycle     = [];
    Trial.Segment(i).T.lcycle     = [];
    Trial.Segment(i).Euler.smooth = [];
    Trial.Segment(i).Euler.rcycle = [];
    Trial.Segment(i).Euler.lcycle = [];
    Trial.Segment(i).Euler.units  = 'rad';
    Trial.Segment(i).sequence     = '';
end