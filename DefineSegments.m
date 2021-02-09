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

function Trial = DefineSegments(Participant,Trial)

% -------------------------------------------------------------------------
% Pelvis parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
RASI = permute(Trial.Marker(1).Trajectory.smooth,[2,3,1]);
RILC = permute(Trial.Marker(2).Trajectory.smooth,[2,3,1]);
RPSI = permute(Trial.Marker(3).Trajectory.smooth,[2,3,1]);
LPSI = permute(Trial.Marker(4).Trajectory.smooth,[2,3,1]);
LILC = permute(Trial.Marker(5).Trajectory.smooth,[2,3,1]);
LASI = permute(Trial.Marker(6).Trajectory.smooth,[2,3,1]);
S1   = permute(Trial.Marker(7).Trajectory.smooth,[2,3,1]);
% Pelvis axes (Dumas and Wojtusch 2018)
Z5       = Vnorm_array3(RASI-LASI);
Y5       = Vnorm_array3(cross(RASI-(RPSI+LPSI)/2, ...
                              LASI-(RPSI+LPSI)/2));
X5       = Vnorm_array3(cross(Y5,Z5));
w_pelvis = mean(sqrt(sum((RASI-LASI).^2))); % pelvis width
% Determination of the lumbar joint centre by regression (Dumas and Wojtusch 2018)
if strcmp(Participant.gender,'Female')
    LJC(1) = -34.0/100;
    LJC(2) = 4.9/100;
    LJC(3) = 0.0/100;
elseif strcmp(Participant.gender,'Male')
    LJC(1) = -33.5/100;
    LJC(2) = -3.2/100;
    LJC(3) = 0.0/100;
end
LJC = (RASI+LASI)/2 + ...
      LJC(1)*w_pelvis*X5 + LJC(2)*w_pelvis*Y5 + LJC(3)*w_pelvis*Z5;
% Store virtual marker
Trial.Vmarker(9).label             = 'LJC';
Trial.Vmarker(9).Trajectory.smooth = permute(LJC,[3,1,2]);  
% Determination of the hip joint centre by regression (Dumas and Wojtusch 2018)
if strcmp(Participant.gender,'Female')
    R_HJC(1) = -13.9/100;
    R_HJC(2) = -33.6/100;
    R_HJC(3) = 37.2/100;
    L_HJC(1) = -13.9/100;
    L_HJC(2) = -33.6/100;
    L_HJC(3) = -37.2/100;
elseif strcmp(Participant.gender,'Male')
    R_HJC(1) = -9.5/100;
    R_HJC(2) = -37.0/100;
    R_HJC(3) = 36.1/100;
    L_HJC(1) = -9.5/100;
    L_HJC(2) = -37.0/100;
    L_HJC(3) = -36.1/100;
end
RHJC = (RASI+LASI)/2 + ...
       R_HJC(1)*w_pelvis*X5 + R_HJC(2)*w_pelvis*Y5 + R_HJC(3)*w_pelvis*Z5;
LHJC = (RASI+LASI)/2 + ...
       L_HJC(1)*w_pelvis*X5 + L_HJC(2)*w_pelvis*Y5 + L_HJC(3)*w_pelvis*Z5;
% Store virtual markers
Trial.Vmarker(4).label             = 'RHJC';
Trial.Vmarker(4).Trajectory.smooth = permute(RHJC,[3,1,2]);  
Trial.Vmarker(8).label             = 'LHJC';
Trial.Vmarker(8).Trajectory.smooth = permute(LHJC,[3,1,2]);  
% Pelvis parameters (Dumas and Ch�ze 2007) = Pelvis duplicated into
% Segments 5 and 9 to manage right and left lower limb kinematic chains
rP5                         = LJC;
rD5                         = (RHJC+LHJC)/2;
w5                          = Z5;
u5                          = X5;
Trial.Segment(5).Q.smooth   = [u5;rP5;rD5;w5];
Trial.Segment(5).rM.smooth  = [RASI,LASI,RPSI,LPSI,RILC,LILC,S1];
Trial.Segment(10).Q.smooth  = [u5;rP5;rD5;w5];
Trial.Segment(10).rM.smooth = [RASI,LASI,RPSI,LPSI,RILC,LILC,S1];

% -------------------------------------------------------------------------
% Right femur parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
RGTR = permute(Trial.Marker(8).Trajectory.smooth,[2,3,1]);
RTHI = permute(Trial.Marker(9).Trajectory.smooth,[2,3,1]);
RKNE = permute(Trial.Marker(10).Trajectory.smooth,[2,3,1]);
RKNM = permute(Trial.Marker(11).Trajectory.smooth,[2,3,1]);
% Knee joint centre
RKJC = (RKNE+RKNM)/2;
% Store virtual marker
Trial.Vmarker(3).label             = 'RKJC';
Trial.Vmarker(3).Trajectory.smooth = permute(RKJC,[3,1,2]);  
% Femur axes (Dumas and Wojtusch 2018)
Y4 = Vnorm_array3(RHJC-RKJC);
X4 = Vnorm_array3(cross(RKNE-RHJC,RKJC-RHJC));
Z4 = Vnorm_array3(cross(X4,Y4));
% Femur parameters (Dumas and Ch�ze 2007)
rP4                        = RHJC;
rD4                        = RKJC;
w4                         = Z4;
u4                         = X4;
Trial.Segment(4).Q.smooth  = [u4;rP4;rD4;w4];
Trial.Segment(4).rM.smooth = [RGTR,RHJC,RKNE,RKNM]; % Not sure about constant orientation of RTHI

% -------------------------------------------------------------------------
% Right Tibia/fibula parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
RFAX = permute(Trial.Marker(12).Trajectory.smooth,[2,3,1]);
RTTA = permute(Trial.Marker(13).Trajectory.smooth,[2,3,1]);
RTIB = permute(Trial.Marker(14).Trajectory.smooth,[2,3,1]);
RANK = permute(Trial.Marker(15).Trajectory.smooth,[2,3,1]);
RMED = permute(Trial.Marker(16).Trajectory.smooth,[2,3,1]);
% Ankle joint centre
RAJC = (RANK+RMED)/2;
% Store virtual marker
Trial.Vmarker(2).label             = 'RAJC';
Trial.Vmarker(2).Trajectory.smooth = permute(RAJC,[3,1,2]);  
% Tibia/fibula axes (Dumas and Wojtusch 2018)
Y3 = Vnorm_array3(RKJC-RAJC);
X3 = Vnorm_array3(cross(RAJC-RFAX,RKJC-RFAX));
Z3 = Vnorm_array3(cross(X3,Y3));
% Tibia/fibula parameters (Dumas and Ch�ze 2007)
rP3                        = RKJC;
rD3                        = RAJC;
w3                         = Z3;
u3                         = X3;
Trial.Segment(3).Q.smooth  = [u3;rP3;rD3;w3];
Trial.Segment(3).rM.smooth = [RFAX,RTTA,RANK,RMED]; % Not sure about constant orientation of RTIB

% -------------------------------------------------------------------------
% Right foot parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
RHEE = permute(Trial.Marker(17).Trajectory.smooth,[2,3,1]);
RFMH = permute(Trial.Marker(19).Trajectory.smooth,[2,3,1]);
RVMH = permute(Trial.Marker(21).Trajectory.smooth,[2,3,1]);
% Metatarsal joint centre (Dumas and Wojtusch 2018)
RMJC = (RFMH+RVMH)/2;
% Store virtual marker
Trial.Vmarker(1).label             = 'RMJC';
Trial.Vmarker(1).Trajectory.smooth = permute(RMJC,[3,1,2]);  
% Foot axes (Dumas and Wojtusch 2018)
X2 = Vnorm_array3(RMJC-RHEE);
Y2 = Vnorm_array3(cross(RVMH-RHEE,RFMH-RHEE));
Z2 = Vnorm_array3(cross(X2,Y2));
% Foot parameters (Dumas and Ch�ze 2007)
rP2                        = RAJC;
rD2                        = RMJC;
w2                         = Z2;
u2                         = X2;
Trial.Segment(2).Q.smooth  = [u2;rP2;rD2;w2];
Trial.Segment(2).rM.smooth = [RHEE,RFMH,RVMH];

% -------------------------------------------------------------------------
% Left femur parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
LGTR = permute(Trial.Marker(23).Trajectory.smooth,[2,3,1]);
LTHI = permute(Trial.Marker(24).Trajectory.smooth,[2,3,1]);
LKNE = permute(Trial.Marker(25).Trajectory.smooth,[2,3,1]);
LKNM = permute(Trial.Marker(26).Trajectory.smooth,[2,3,1]);
% Knee joint centre
LKJC = (LKNE+LKNM)/2;
% Store virtual marker
Trial.Vmarker(7).label             = 'LKJC';
Trial.Vmarker(7).Trajectory.smooth = permute(LKJC,[3,1,2]);  
% Femur axes (Dumas and Wojtusch 2018)
Y9 = Vnorm_array3(LHJC-LKJC);
X9 = -Vnorm_array3(cross(LKNE-LHJC,LKJC-LHJC));
Z9 = Vnorm_array3(cross(X9,Y9));
% Femur parameters (Dumas and Ch�ze 2007)
rP9                        = LHJC;
rD9                        = LKJC;
w9                         = Z9;
u9                         = X9;
Trial.Segment(9).Q.smooth  = [u9;rP9;rD9;w9];
Trial.Segment(9).rM.smooth = [LGTR,LHJC,LKNE,LKNM]; % Not sure about constant orientation of LTHI

% -------------------------------------------------------------------------
% Left Tibia/fibula parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
LFAX = permute(Trial.Marker(27).Trajectory.smooth,[2,3,1]);
LTTA = permute(Trial.Marker(28).Trajectory.smooth,[2,3,1]);
LTIB = permute(Trial.Marker(29).Trajectory.smooth,[2,3,1]);
LANK = permute(Trial.Marker(30).Trajectory.smooth,[2,3,1]);
LMED = permute(Trial.Marker(31).Trajectory.smooth,[2,3,1]);
% Ankle joint centre
LAJC = (LANK+LMED)/2;
% Store virtual marker
Trial.Vmarker(6).label             = 'LAJC';
Trial.Vmarker(6).Trajectory.smooth = permute(LAJC,[3,1,2]);  
% Tibia/fibula axes (Dumas and Wojtusch 2018)
Y8 = Vnorm_array3(LKJC-LAJC);
X8 = -Vnorm_array3(cross(LAJC-LFAX,LKJC-LFAX));
Z8 = Vnorm_array3(cross(X8,Y8));
% Tibia/fibula parameters (Dumas and Ch�ze 2007)
rP8                        = LKJC;
rD8                        = LAJC;
w8                         = Z8;
u8                         = X8;
Trial.Segment(8).Q.smooth  = [u8;rP8;rD8;w8];
Trial.Segment(8).rM.smooth = [LFAX,LTTA,LANK,LMED]; % Not sure about constant orientation of LTIB

% -------------------------------------------------------------------------
% Left foot parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
LHEE = permute(Trial.Marker(32).Trajectory.smooth,[2,3,1]);
LFMH = permute(Trial.Marker(34).Trajectory.smooth,[2,3,1]);
LVMH = permute(Trial.Marker(36).Trajectory.smooth,[2,3,1]);
% Metatarsal joint centre (Dumas and Wojtusch 2018)
LMJC = (LFMH+LVMH)/2;
% Store virtual marker
Trial.Vmarker(5).label             = 'LMJC';
Trial.Vmarker(5).Trajectory.smooth = permute(LMJC,[3,1,2]);  
% Foot axes (Dumas and Wojtusch 2018)
X7 = Vnorm_array3(LMJC-LHEE);
Y7 = -Vnorm_array3(cross(LVMH-LHEE,LFMH-LHEE));
Z7 = Vnorm_array3(cross(X7,Y7));
% Foot parameters (Dumas and Ch�ze 2007)
rP7                        = LAJC;
rD7                        = LMJC;
w7                         = Z7;
u7                         = X7;
Trial.Segment(7).Q.smooth  = [u7;rP7;rD7;w7];
Trial.Segment(7).rM.smooth = [LHEE,LFMH,LVMH];

% -------------------------------------------------------------------------
% Lumbar parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
CLAV = permute(Trial.Marker(48).Trajectory.smooth,[2,3,1]);
STRN = permute(Trial.Marker(49).Trajectory.smooth,[2,3,1]);
C7   = permute(Trial.Marker(41).Trajectory.smooth,[2,3,1]);
T8   = permute(Trial.Marker(45).Trajectory.smooth,[2,3,1]);
T12  = permute(Trial.Marker(47).Trajectory.smooth,[2,3,1]);
L1   = permute(Trial.Marker(38).Trajectory.smooth,[2,3,1]);
% Thorax width (Dumas and Wojtusch 2018)
W11 = mean(sqrt(sum((STRN-T8).^2)));
% Determination of thoracic joint centre by regression (Dumas and Wojtusch 2018)
tX11 = Vnorm_array3(STRN-T8);
tY11 = Vnorm_array3(T8-T12);
tZ11 = cross(tX11,tY11);
tX11 = cross(tY11,tZ11);
if strcmp(Participant.gender,'Female')
    angle = 92;
    coeff = 0.50;
elseif strcmp(Participant.gender,'Male')
    angle = 94;
    coeff = 0.52;
end
R11 = [cosd(angle) sind(angle) 0 0; ...
       -sind(angle) cosd(angle) 0 0;
       0 0 1 0; ...
       0 0 0 1];
TJC = Mprod_array3(Mprod_array3([tX11 tY11 tZ11 T12; ...
                   repmat([0 0 0 1],[1,1,size(T12,3)])], ...
                   repmat(R11,[1,1,size(T12,3)])), ...
                   repmat([0; coeff*W11; 0; 1],[1,1,size(T12,3)]));
TJC = TJC(1:3,:,:);
% Store virtual marker
Trial.Vmarker(11).label             = 'TJC';
Trial.Vmarker(11).Trajectory.smooth = permute(TJC,[3,1,2]);  
% Lumbar axes (Dumas and Wojtusch 2018)
Y11 = Vnorm_array3(TJC-LJC);
Z11 = Z5; % no axial rotation at lumbar joint centre assumed
X11 = Vnorm_array3(cross(Y11,Z11));
% Lumbar parameters
rP11                        = TJC;
rD11                        = LJC;
w11                         = Z11;
u11                         = X11;
Trial.Segment(11).Q.smooth  = [u11;rP11;rD11;w11];
Trial.Segment(11).rM.smooth = [TJC,RPSI,LPSI]; % no axial rotation at lumbar joint centre assumed

% -------------------------------------------------------------------------
% Thorax parameters
% -------------------------------------------------------------------------
% Thorax width (Dumas and Wojtusch 2018)
W12 = mean(sqrt(sum((CLAV-C7).^2)));
% Determination of the cervical joint centre by regression (Dumas and Wojtusch 2018)
tX12 = Vnorm_array3(CLAV-C7);
tZ12 = Vnorm_array3(cross(STRN-C7,CLAV-C7));
tY12 = Vnorm_array3(cross(tZ12,tX12));
if strcmp(Participant.gender,'Female')
    angle = -14;
    coeff = 0.53;
elseif strcmp(Participant.gender,'Male')
    angle = -8;
    coeff = 0.55;
end
R12 = [cosd(angle) sind(angle) 0 0; ...
       -sind(angle) cosd(angle) 0 0;
       0 0 1 0; ...
       0 0 0 1];
CJC = Mprod_array3(Mprod_array3([tX12 tY12 tZ12 C7; ...
                   repmat([0 0 0 1],[1,1,size(C7,3)])], ...
                   repmat(R12,[1,1,size(C7,3)])), ...
                   repmat([coeff*W12; 0; 0; 1],[1,1,size(C7,3)]));
CJC = CJC(1:3,:,:);
% Store virtual marker
Trial.Vmarker(12).label             = 'CJC';
Trial.Vmarker(12).Trajectory.smooth = permute(CJC,[3,1,2]);  
% Thorax axes (Dumas and Ch�ze 2007)
Y12 = Vnorm_array3(CJC-TJC);
Z12 = Vnorm_array3(cross(CLAV-TJC,CJC-TJC));
X12 = Vnorm_array3(cross(Y12,Z12));
% Thorax parameters
rP12                        = CJC;
rD12                        = TJC;
w12                         = Z12;
u12                         = X12;
Trial.Segment(12).Q.smooth  = [u12;rP12;rD12;w12];
Trial.Segment(12).rM.smooth = [CLAV,C7,STRN,TJC];

% -------------------------------------------------------------------------
% Head with neck parameters
% -------------------------------------------------------------------------
% Extract marker trajectories
RFHD = permute(Trial.Marker(50).Trajectory.smooth,[2,3,1]);
RBHD = permute(Trial.Marker(51).Trajectory.smooth,[2,3,1]);
LFHD = permute(Trial.Marker(52).Trajectory.smooth,[2,3,1]);
LBHD = permute(Trial.Marker(53).Trajectory.smooth,[2,3,1]);
% Head vertex (Dumas and Wojtusch 2018)
VER = (RFHD+RBHD+LFHD+LBHD)/4; % assimilated to the head vertex described in Dumas and Wojtusch 2018
% Store virtual marker
Trial.Vmarker(12).label             = 'VER';
Trial.Vmarker(12).Trajectory.smooth = permute(VER,[3,1,2]);  
% Head axes
Y13 = Vnorm_array3(VER-CJC);
Z13 = Vnorm_array3(cross((RFHD+LFHD)/2-CJC,(RBHD+LBHD)/2-CJC));
X13 = Vnorm_array3(cross(Y13,Z13));
% Head parameters
rP13                        = VER;
rD13                        = CJC;
w13                         = Z13;
u13                         = X13;
Trial.Segment(13).Q.smooth  = [u13;rP13;rD13;w13];
Trial.Segment(13).rM.smooth = [RFHD,RBHD,LFHD,LBHD];