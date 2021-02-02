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

function Trial = InitialiseEMGSignals(Trial,EMG)

% Set EMGset
EMGSet = {'R_RA','L_RA',...
          'R_EO','L_EO',...
          'R_LES','L_LES',...
          'R_ESI','L_ESI',...
          'R_MTF','L_MTF',...
          'R_GMED','L_GMED',...
          'R_RF','L_RF',...
          'R_SM','L_SM'};      
         
% Set related task
trialList = {'Endurance_Ito','Endurance_Ito',...
             'Endurance_Ito','Endurance_Ito',...
             'Endurance_Sorensen','Endurance_Sorensen',...
             'Endurance_Sorensen','Endurance_Sorensen',...
             'Endurance_Sorensen','Endurance_Sorensen',...
             'sMVC_R_Gmed','sMVC_L_Gmed',...
             'sMVC_R_Rf','sMVC_L_Rf',...
             'sMVC_R_Semiten','sMVC_L_Semiten'};   

% Initialise EMGs
for i = 1:length(EMGSet)
    Trial.EMG(i).label              = EMGSet{i};
    if isfield(EMG,EMGSet{i})
        Trial.EMG(i).Signal.raw     = EMG.(EMGSet{i});
        Trial.EMG(i).Signal.filt    = [];
        Trial.EMG(i).Signal.rect    = [];
        Trial.EMG(i).Signal.smooth  = [];
        Trial.EMG(i).Signal.norm    = [];
        Trial.EMG(i).Signal.cycle   = []; %smooth
        Trial.EMG(i).Signal.rcycle  = []; %smooth
        Trial.EMG(i).Signal.lcycle  = []; %smooth
        Trial.EMG(i).Signal.cyclen  = []; %norm
        Trial.EMG(i).Signal.rcyclen = []; %norm
        Trial.EMG(i).Signal.lcyclen = []; %norm
    else
        Trial.EMG(i).Signal.raw     = [];
        Trial.EMG(i).Signal.filt    = [];
        Trial.EMG(i).Signal.rect    = [];
        Trial.EMG(i).Signal.smooth  = [];
        Trial.EMG(i).Signal.norm    = [];
        Trial.EMG(i).Signal.cycle   = []; %smooth
        Trial.EMG(i).Signal.rcycle  = []; %smooth
        Trial.EMG(i).Signal.lcycle  = []; %smooth
        Trial.EMG(i).Signal.cyclen  = []; %norm
        Trial.EMG(i).Signal.rcyclen = []; %norm
        Trial.EMG(i).Signal.lcyclen = []; %norm
    end
    Trial.EMG(i).Processing.filt       = 'none';
    Trial.EMG(i).Processing.smooth     = 'none';
    Trial.EMG(i).Processing.normMethod = 'none';
    Trial.EMG(i).Processing.normTask   = trialList{i};
    Trial.EMG(i).Processing.normValue  = [];
end