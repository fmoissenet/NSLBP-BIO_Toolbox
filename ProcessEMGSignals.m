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
% Dependencies : - Biomechanical Toolkit (BTK): https://github.com/Biomechanical-ToolKit/BTKCore
% -------------------------------------------------------------------------
% This work is licensed under the Creative Commons Attribution - 
% NonCommercial 4.0 International License. To view a copy of this license, 
% visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% -------------------------------------------------------------------------

function [Calibration,Trial] = ProcessEMGSignals(Calibration,Trial,defCal,fmethod,smethod,nmethod)

for i = 1:size(Trial.EMG,2)
    if ~isempty(Trial.EMG(i).Signal.raw)
        
        % -----------------------------------------------------------------
        % REPLACE NANs BY ZEROs
        % -----------------------------------------------------------------   
        Trial.EMG(i).Signal.filt                                  = Trial.EMG(i).Signal.raw;
        Trial.EMG(i).Signal.filt(isnan(Trial.EMG(i).Signal.filt)) = 0;

        % -----------------------------------------------------------------
        % ZEROING AND FILTER EMG SIGNAL
        % -----------------------------------------------------------------
        Trial.EMG(i).Signal.filt = Trial.EMG(i).Signal.filt-nanmean(Trial.EMG(i).Signal.filt);
        
        % Method 1: No filtering
        if strcmp(fmethod.type,'none')
            Trial.EMG(i).Signal.filt     = Trial.EMG(i).Signal.filt;
            Trial.EMG(i).Processing.filt = fmethod.type;
        
        % Method 2: Band pass filter (Butterworth 4nd order, [fmethod.parameter fmethod.parameter] Hz)
        elseif strcmp(fmethod.type,'butterBand4')
            [B,A]                        = butter(2,[fmethod.parameter(1) fmethod.parameter(2)]./(Trial.fanalog/2),'bandpass');
            Trial.EMG(i).Signal.filt     = filtfilt(B,A,Trial.EMG(i).Signal.filt);
            Trial.EMG(i).Processing.filt = fmethod.type;
        end

        % -----------------------------------------------------------------
        % RECTIFY EMG SIGNAL
        % -----------------------------------------------------------------
        Trial.EMG(i).Signal.rect = abs(Trial.EMG(i).Signal.filt);
        
        % -----------------------------------------------------------------
        % SMOOTH EMG SIGNAL
        % -----------------------------------------------------------------        
        % Method 1: No smoothing
        if strcmp(smethod.type,'none')
            Trial.EMG(i).Signal.smooth     = Trial.EMG(i).Signal.rect;
            Trial.EMG(i).Processing.smooth = smethod.type;
        
        % Method 2: Low pass filter (Butterworth 2nd order, [smethod.parameter] Hz)
        elseif strcmp(smethod.type,'butterLow2')
            [B,A]                          = butter(1,smethod.parameter/(Trial.fanalog/2),'low');
            Trial.EMG(i).Signal.smooth     = filtfilt(B,A,Trial.EMG(i).Signal.rect);
            Trial.EMG(i).Processing.smooth = smethod.type;
        
        % Method 3: Moving average (window of [smethod.parameter] frames)
        elseif strcmp(smethod.type,'movmean')
            Trial.EMG(i).Signal.smooth = smoothdata(Trial.EMG(i).Signal.rect,'movmean',smethod.parameter);
            Trial.EMG(i).Processing.smooth = 'movmean';
        
        % Method 4: Moving average (window of [smethod.parameter] frames)
        elseif strcmp(smethod.type,'movmedian')
            Trial.EMG(i).Signal.smooth = smoothdata(Trial.EMG(i).Signal.rect,'movmedian',smethod.parameter);
            Trial.EMG(i).Processing.smooth = 'movmedian';
        
        % Method 5: Signal root mean square (RMS) (window of [smethod.parameter] frames)
        elseif strcmp(smethod.type,'rms')
            Trial.EMG(i).Signal.smooth     = envelope(Trial.EMG(i).Signal.filt,smethod.parameter,'rms')
            Trial.EMG(i).Processing.smooth = 'rms';
        end
        
        % -----------------------------------------------------------------
        % APPLY SIGNAL CALIBRATION
        % -----------------------------------------------------------------        
        if defCal ~= 1
            if ~isempty(Calibration.EMG(i).normValue)
                Trial.EMG(i).Signal.norm     = Trial.EMG(i).Signal.smooth/Calibration.EMG(i).normValue;
                Trial.EMG(i).Processing.norm = Calibration.EMG(i).normMethod;
            end
        end

        % -----------------------------------------------------------------
        % KEEP ONLY MEAN VALUES FOR NORMALISATION DATA
        % -----------------------------------------------------------------      
        if defCal == 1
            
            % Method 1: No normalisation
            if strcmp(nmethod.type,'none')
                Calibration.EMG(i).normValue  = 1;
                Calibration.EMG(i).normMethod = nmethod.type;
            
            % Method 2: Mean value during subMVC task
            elseif strcmp(nmethod.type,'sMVC')
                if ~isempty(strfind(Trial.file,Trial.EMG(i).Processing.normTask))                    
                    start                         = Trial.Event(1).value*Trial.fanalog/Trial.fmarker;
                    stop                          = Trial.Event(2).value*Trial.fanalog/Trial.fmarker;
                    Calibration.EMG(i).label      = Trial.EMG(i).label;
                    Calibration.EMG(i).file       = Trial.file;
                    temp = [];
                    for j = 1:size(start,2)
                        temp = [temp; Trial.EMG(i).Signal.smooth(start(j):stop(j),1)];
                    end
                    Calibration.EMG(i).normValue  = mean(temp);
                    Calibration.EMG(i).normMethod = nmethod.type;
                    clear temp;
                end                
            end
            
        end        
    end
end