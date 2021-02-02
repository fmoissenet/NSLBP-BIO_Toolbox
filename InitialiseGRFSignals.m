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

function Trial = InitialiseGRFSignals(Trial,GRF,GRFmeta)

if ~isempty(GRF)
    for i = 1:size(GRF,1)
        
        % Store plateforme corners location
        Trial.GRF(i).label         = ['Forceplate_',num2str(i)];
        Trial.GRF(i).side          = '';
        Trial.GRF(i).cycle         = [];
        Trial.GRF(i).Location      = [];
        corners                    = GRFmeta.children.CORNERS.info.values(:,:,i);
        Trial.GRF(i).Location.X(1) = min([corners(1,1) corners(1,2) corners(1,3) corners(1,4)])*1e-3; % Convert mm to m
        Trial.GRF(i).Location.X(2) = max([corners(1,1) corners(1,2) corners(1,3) corners(1,4)])*1e-3; % Convert mm to m
        Trial.GRF(i).Location.Y(1) = min([corners(2,1) corners(2,2) corners(2,3) corners(2,4)])*1e-3; % Convert mm to m
        Trial.GRF(i).Location.Y(2) = max([corners(2,1) corners(2,2) corners(2,3) corners(2,4)])*1e-3; % Convert mm to m
        clear corners;
                
        % Initialise CoP, force and moment
        Trial.GRF(i).Signal.P.raw      = [];
        Trial.GRF(i).Signal.P.filt     = [];
        Trial.GRF(i).Signal.P.smooth   = [];
        Trial.GRF(i).Signal.P.cycle    = [];
        Trial.GRF(i).Signal.P.rcycle   = [];
        Trial.GRF(i).Signal.P.lcycle   = [];
        Trial.GRF(i).Signal.F.raw      = [];
        Trial.GRF(i).Signal.F.filt     = [];
        Trial.GRF(i).Signal.F.smooth   = [];
        Trial.GRF(i).Signal.F.cycle    = [];
        Trial.GRF(i).Signal.F.rcycle   = [];
        Trial.GRF(i).Signal.F.lcycle   = [];
        Trial.GRF(i).Signal.M.raw      = [];
        Trial.GRF(i).Signal.M.filt     = [];
        Trial.GRF(i).Signal.M.smooth   = [];
        Trial.GRF(i).Signal.M.cycle    = [];
        Trial.GRF(i).Signal.M.rcycle   = [];
        Trial.GRF(i).Signal.M.lcycle   = [];
        Trial.GRF(i).Processing.filt   = 'none';
        Trial.GRF(i).Processing.smooth = 'none';  
    end
end