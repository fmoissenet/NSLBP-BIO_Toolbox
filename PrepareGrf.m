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
% Description  : This routine aims to load and process ground reaction
%                forces and moments as well as centre of pressure.
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

TODO


function Grf = PrepareGrf(Grf,n0,f)

if noGrf == 0
    % Modify the ICS
    if strcmp(system,'Qualisys')
        if strcmp(side,'Right') 
            Grf.P = [Grf.P(:,1) Grf.P(:,3) -Grf.P(:,2)];
            Grf.F = [Grf.F(:,1) Grf.F(:,3) -Grf.F(:,2)];
            Grf.M = [zeros(size(Grf.M(:,3))) -Grf.M(:,3) zeros(size(Grf.M(:,3)))];
        elseif strcmp(side,'Left') 
            Grf.P = [Grf.P(:,1) Grf.P(:,3) Grf.P(:,2)];
            Grf.F = [Grf.F(:,1) Grf.F(:,3) Grf.F(:,2)];
            Grf.M = [zeros(size(Grf.M(:,3))) Grf.M(:,3) zeros(size(Grf.M(:,3)))];
        end
    end
    % Convert from Nmm to Nm
    if strcmp(btkGetPointsUnit(file,'marker'),'mm')
        Grf.M = Grf.M*10^(-3);
        Grf.P = Grf.P*10^(-3);
    end
    % Filt data (15Hz, 2th order Btw)
    [B,A] = butter(2,15/(f2/2),'low');
    Grf.F(:,1) = filtfilt(B,A,Grf.F(:,1));
    Grf.F(:,2) = filtfilt(B,A,Grf.F(:,2));
    Grf.F(:,3) = filtfilt(B,A,Grf.F(:,3));
    Grf.M(:,1) = filtfilt(B,A,Grf.M(:,1));
    Grf.M(:,2) = filtfilt(B,A,Grf.M(:,2));
    Grf.M(:,3) = filtfilt(B,A,Grf.M(:,3));
    % Apply a 10N threshold
    threshold = 10;
    for j = 1:n2
        if Grf.F(j,2) < threshold;
            Grf.P(j,:) = zeros(1,3);
            Grf.F(j,:) = zeros(1,3);
            Grf.M(j,:) = zeros(1,3);
        end
    end
    % Resample to the marker trajectories number of frames
    x = 1:n2;
    xx = linspace(1,n2,n1);
    if ~isnan(mean(Grf.F(1,:,:)))
        Grf.P = interp1(x,Grf.P,xx,'spline');
        Grf.F = interp1(x,Grf.F,xx,'spline');
        Grf.M = interp1(x,Grf.M,xx,'spline');
    end
    % Set walking direction
    Grf.P(:,1) = minusX*Grf.P(:,1);
    Grf.P(:,3) = minusX*Grf.P(:,3);
    Grf.F(:,1) = minusX*Grf.F(:,1);
    Grf.F(:,3) = minusX*Grf.F(:,3);
    Grf.M(:,1) = minusX*Grf.M(:,1);
    Grf.M(:,3) = minusX*Grf.M(:,3);
    % Store Marker as 3-array vectors
    Grf.P = permute(Grf.P,[2,3,1]);
    Grf.F = permute(Grf.F,[2,3,1]);
    Grf.M = permute(Grf.M,[2,3,1]);
end