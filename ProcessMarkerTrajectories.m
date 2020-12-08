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
% Description  : This routine aims to process 3D marker trajectories
% Inputs       : To be defined
% Outputs      : To be defined
% -------------------------------------------------------------------------
% Dependencies : - Biomechanical Toolkit (BTK): https://github.com/Biomechanical-ToolKit/BTKCore
%                - PredictMissingMarkers: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0152616
% -------------------------------------------------------------------------
% This work is licensed under the Creative Commons Attribution - 
% NonCommercial 4.0 International License. To view a copy of this license, 
% visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% -------------------------------------------------------------------------

function Trial = ProcessMarkerTrajectories(Static,Trial,fmethod,smethod)

% -------------------------------------------------------------------------
% IDENTIFY MISSING TRAJECTORIES AND GAPS IN TRAJECTORIES 
% -------------------------------------------------------------------------
for i = 1:size(Trial.Marker,2)

    % Special case of static
    if isempty(Static)
        if isempty(Trial.Marker(i).Trajectory.raw)
            Trial.Marker(i).Trajectory.fill = [];
        else
            Trial.Marker(i).Trajectory.fill = mean(Trial.Marker(i).Trajectory.raw,1);
        end
        Trial.n0                        = 1;
        Trial.n1                        = 1;  
    end
    
    % Missing marker trajectory
    if isempty(Trial.Marker(i).Trajectory.raw)
        Trial.Marker(i).Gap(1).frames         = 1:Trial.n1;
        Trial.Marker(i).Gap(1).reconstruction = 'none';
        Trial.Marker(i).smoothing             = 'none';

    % Marker trajectory with gaps
    elseif ~isempty(Trial.Marker(i).Trajectory.raw)

        % Replace [0 0 0] by NaN
        for j = 1:Trial.n1
            if Trial.Marker(i).Trajectory.raw(j,:) == [0 0 0]
               Trial.Marker(i).Trajectory.fill(j,:) = nan(1,3);
            else
               Trial.Marker(i).Trajectory.fill(j,:) = Trial.Marker(i).Trajectory.raw(j,:);
            end
        end

        % Find gaps
        start = 0;
        stop  = 0;
        k     = 0;
        for j = 1:Trial.n1-1
            if isnan(Trial.Marker(i).Trajectory.fill(j,:))
                if start == 0
                    start = j;
                end
                if ~isnan(Trial.Marker(i).Trajectory.fill(j+1,:))
                    if start ~= 0
                        stop   = j;
                        k      = k+1;    
                        Trial.Marker(i).Gap(k).frames         = start:stop;
                        Trial.Marker(i).Gap(k).reconstruction = 'none';
                        start  = 0;
                        stop   = 0;
                    end
                elseif j+1 == Trial.n1
                    if isnan(Trial.Marker(i).Trajectory.fill(j+1,:))
                        if start ~= 0
                            stop   = j+1;
                            k      = k+1;    
                            Trial.Marker(i).Gap(k).frames         = start:stop;
                            Trial.Marker(i).Gap(k).reconstruction = 'none';
                            start  = 0;
                            stop   = 0;
                        end   
                    end
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% TRAJECTORIES GAP FILLING (NOT ALLOWED FOR STATIC)
% -------------------------------------------------------------------------
if ~isempty(Static)
    
    % Method 0: None
    if strcmp(fmethod.type,'none')
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                if ~isempty(Trial.Marker(i).Gap)
                    for j = 1:size(Trial.Marker(i).Gap,2)
                        Trial.Marker(i).Gap(j).reconstruction = 'none';
                    end
                end
            end
        end
    end
    
    % Method 1: Linear interpolation
    %           - At least 1 point before and 1 point after gap is required (1/1
    %             are used here)
    if strcmp(fmethod.type,'linear')
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                if ~isempty(Trial.Marker(i).Gap)
                    for j = 1:size(Trial.Marker(i).Gap,2)
                        if size(Trial.Marker(i).Gap(j).frames,2) < fmethod.gapThreshold
                            if Trial.Marker(i).Gap(j).frames(1) > 1 && ...
                               Trial.Marker(i).Gap(j).frames(end) < Trial.n1
                                Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-1: ...
                                                                Trial.Marker(i).Gap(j).frames(end)+1,:) = ...
                                fillmissing(Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-1: ...
                                                                            Trial.Marker(i).Gap(j).frames(end)+1,:),'linear');
                            end
                            Trial.Marker(i).Gap(j).reconstruction = 'linear';
                        end
                    end
                end
            end
        end
    end

    % Method 2: Cubic spline interpolation
    %           - At least 2 point before and 2 point after gap is required 
    %             (10/10 are used here)
    if strcmp(fmethod.type,'spline')
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                if ~isempty(Trial.Marker(i).Gap)
                    for j = 1:size(Trial.Marker(i).Gap,2)
                        if size(Trial.Marker(i).Gap(j).frames,2) < fmethod.gapThreshold
                            if Trial.Marker(i).Gap(j).frames(1) > 10 && ...
                               Trial.Marker(i).Gap(j).frames(end) < Trial.n1-9
                                Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-10: ...
                                                                Trial.Marker(i).Gap(j).frames(end)+10,:) = ...
                                fillmissing(Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-10: ...
                                                                            Trial.Marker(i).Gap(j).frames(end)+10,:),'spline');
                            end
                            Trial.Marker(i).Gap(j).reconstruction = 'spline';
                        end
                    end
                end
            end
        end
    end

    % Method 3: Shape-preserving piecewise cubic interpolation
    %           - At least 2 point before and 2 point after gap is required 
    %             (10/10 are used here)
    if strcmp(fmethod.type,'pchip')
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                if ~isempty(Trial.Marker(i).Gap)
                    for j = 1:size(Trial.Marker(i).Gap,2)
                        if size(Trial.Marker(i).Gap(j).frames,2) < fmethod.gapThreshold
                            if Trial.Marker(i).Gap(j).frames(1) > 10 && ...
                               Trial.Marker(i).Gap(j).frames(end) < Trial.n1-9
                                Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-10: ...
                                                                Trial.Marker(i).Gap(j).frames(end)+10,:) = ...
                                fillmissing(Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-10: ...
                                                                            Trial.Marker(i).Gap(j).frames(end)+10,:),'pchip');
                            end
                            Trial.Marker(i).Gap(j).reconstruction = 'pchip';
                        end
                    end
                end
            end
        end
    end

    % Method 4: Modified Akima cubic Hermite interpolation
    %           - At least 1 point before and 1 point after gap is required 
    %             (10/10 are used here)
    if strcmp(fmethod.type,'makima')
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                if ~isempty(Trial.Marker(i).Gap)
                    for j = 1:size(Trial.Marker(i).Gap,2)
                        if size(Trial.Marker(i).Gap(j).frames,2) < fmethod.gapThreshold
                            if Trial.Marker(i).Gap(j).frames(1) > 10 && ...
                               Trial.Marker(i).Gap(j).frames(end) < Trial.n1-9
                                Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-10: ...
                                                                Trial.Marker(i).Gap(j).frames(end)+10,:) = ...
                                fillmissing(Trial.Marker(i).Trajectory.fill(Trial.Marker(i).Gap(j).frames(1)-10: ...
                                                                            Trial.Marker(i).Gap(j).frames(end)+10,:),'makima');
                            end
                            Trial.Marker(i).Gap(j).reconstruction = 'makima';
                        end
                    end
                end
            end
        end
    end

    % Method 5: Marker trajectories intercorrelation (https://doi.org/10.1371/journal.pone.0152616)
    if strcmp(fmethod.type,'intercor')
        tMarker = [];
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                tMarker = [tMarker Trial.Marker(i).Trajectory.fill];
            end
        end
        tMarker = PredictMissingMarkers(tMarker,'Algorithm',2);
        k = 0;
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                k = k+1;
                Trial.Marker(i).Trajectory.fill = tMarker(:,(3*k)-2:3*k);
            end
        end
        clear k tMarker;
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                if ~isempty(Trial.Marker(i).Gap)
                    for j = 1:size(Trial.Marker(i).Gap,2)
                        Trial.Marker(i).Gap(j).reconstruction = 'intercor';
                    end
                end
            end
        end
    end

    % Method 6: Apply rigid body transformation of the related segment on
    %           missing trajectories
    %           - The missing trajectories must be part of a marker related to a
    %             rigid body
    %           - At least 3 other markers, without gap, are needed on each segment
    if strcmp(fmethod.type,'rigid')
        for i = 1:size(Trial.Marker,2)        
            if ~isempty(Trial.Marker(i).Gap)
                for j = 1:size(Trial.Marker(i).Gap,2)

                    % Markers related to a rigid body
                    if strcmp(Trial.Marker(i).type,'landmark') || ...
                            strcmp(Trial.Marker(i).type,'hybrid-landmark') || ...
                            strcmp(Trial.Marker(i).type,'technical')

                        % Identify all available markers of the same segment
                        % without gap during all frames of the processed gap
                        nsegment = Trial.Marker(i).Body.Segment.label;
                        kmarker = [];
                        if strcmp(nsegment,'none') == 0 % Only for available segments
                            for k = 1:size(Trial.Marker,2)
                                if k ~= i
                                    if strcmp(Trial.Marker(k).Body.Segment.label,nsegment) == 1
                                        if ~isempty(Trial.Marker(k).Trajectory.raw)
                                            if isempty(find(isnan(Trial.Marker(k).Trajectory.fill(Trial.Marker(i).Gap(j).frames,1))))
                                                kmarker = [kmarker k];
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        % If at least 3 markers of the same segment are
                        % available, reconstruct the missing marker
                        if size(kmarker,2) >= 3
                            X = [];
                            for k = 1:size(kmarker,2)
                                X = [X; Static.Marker(kmarker(k)).Trajectory.fill];
                            end
                            for t = Trial.Marker(i).Gap(j).frames
                                Y = [];
                                for k = 1:size(kmarker,2)
                                    Y = [Y; Trial.Marker(kmarker(k)).Trajectory.fill(t,:)];
                                end
                                [R,d,rms] = soder(X,Y);
                                Trial.Marker(i).Trajectory.fill(t,:) = ...
                                    Static.Marker(i).Trajectory.fill*R'+d';
                                clear R d;
                            end
                        end
                        clear segment;
                    end
                    Trial.Marker(i).Gap(j).reconstruction = 'rigid';
                end
            end
        end
    end

    % Method 7: Low dimensional Kalman smoothing (http://dx.doi.org/10.1016/j.jbiomech.2016.04.016)
    %           - A set of frames in which all markers are present is required
    if strcmp(fmethod.type,'kalman')
        tMarker = [];
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                tMarker = [tMarker Trial.Marker(i).Trajectory.fill];
            end
        end
        % Parameters
        tol = 1e-3; % Default: 5e-6;
        sigmaR = 5e-2; % Confidence in the model, default: 5e-2
        % Get nice training data - all sets with full markers
        [r,~] = find(isnan(tMarker));
        bins = setdiff(1:length(tMarker),unique(r));
        Train = tMarker(bins,:);
        % Project Training data into pca space
        %[U,V,l] = pca(Train(:,:));
        [U,S,V] = svd(bsxfun(@minus,Train,mean(Train)));
        mPCA = mean(Train(:,:));    
        l = diag(S);
        % Determine required number dimensions for model fitting
        d = find(abs(cumsum(l)./sum(l)-1) < tol,1,'first');
        % Find model process noise
        sigma_a = var(diff(Train));  % Confidence in the measurements
        Q = V(:,1:d)'*((diag(sigma_a)).^2)*V(:,1:d);
        % Forward stage
        Estimate = zeros(length(tMarker),size(tMarker,2));
        frate = 0;
        bin = 1;
        state_pred{1} = randn(d,1);
        state{1} = randn(d,1);
        cov{1} = 1e12*eye(d);
        cov_pred{1} = 1e12*eye(d);
        for j = 2:length(tMarker)+1
            tic
            % Construct measurement matrix
            H = diag(~isnan(tMarker(j-1,:)));
            H(sum(H,2)==0,:) = [];
            % Measurement noise
            R = sigmaR*eye(size(H,1));
            Ht = H*V(:,1:d);        
            % Extract valid measurements
            z = tMarker(j-1,~isnan(tMarker(j-1,:)))';        
            state_pred{j} = state{j-1};
            cov_pred{j} = cov{j-1} + Q;        
            K = cov_pred{j}*Ht'*inv(Ht*cov_pred{j}*Ht' + R);        
            state{j} = state_pred{j} + K*(z - (Ht*state_pred{j} + H*mPCA'));        
            cov{j} = (eye(d) - K*Ht)*cov_pred{j};        
            est = V(:,1:d)*state{j} + mPCA';        
            frate = frate+toc;        
            bin = bin + 1;
        end
        % Backward stage
        state_new = cell(length(state)-1,1);
        cov_new = cell(length(cov)-1,1);
        state_new{end} = state{end};
        cov_new{end} = cov{end};
        frate = 0;
        for j = length(state)-1:-1:2
            tic        
            state_new{j} = state{j} + cov{j}*inv(cov_pred{j})*(state{j+1} - state_pred{j+1});
            cov_new{j} = cov{j} + cov{j}*inv(cov_pred{j})*(cov{j+1} - cov_pred{j+1})*cov{j};        
            Estimate(j-1,:) = V(:,1:d)*state_new{j} + mPCA';        
            frate = frate + toc;          
        end
        Estimate(end,:) = V(:,1:d)*state_new{end} + mPCA';
        k = 0;
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                k = k+1;
                Trial.Marker(i).Trajectory.fill = Estimate(:,(3*k)-2:3*k);
            end
        end
        clear k tMarker;  
        for i = 1:size(Trial.Marker,2)
            if ~isempty(Trial.Marker(i).Trajectory.raw)
                if ~isempty(Trial.Marker(i).Gap)
                    for j = 1:size(Trial.Marker(i).Gap,2)
                        Trial.Marker(i).Gap(j).reconstruction = 'kalman';
                    end
                end
            end
        end
    end
end

% -------------------------------------------------------------------------
% MISSING TRAJECTORIES RECONSTRUCTION (NOT ALLOWED FOR STATIC)
% -------------------------------------------------------------------------
if ~isempty(Static)
    for i = 1:size(Trial.Marker,2)
        if isempty(Trial.Marker(i).Trajectory.raw)
            
            % Markers related to a rigid body (landmarks and hybrid-landmarks)
            if strcmp(Trial.Marker(i).type,'landmark') || ...
               strcmp(Trial.Marker(i).type,'hybrid-landmark') || ...
               strcmp(Trial.Marker(i).type,'technical')
                
                % Identify all available markers of the same segment
                % without gap during all frames of the processed gap
                nsegment = Trial.Marker(i).Body.Segment.label;
                kmarker = [];
                if strcmp(nsegment,'none') == 0 % Only for available segments
                    for k = 1:size(Trial.Marker,2)
                        if k ~= i
                            if strcmp(Trial.Marker(k).Body.Segment.label,nsegment) == 1
                                if ~isempty(Trial.Marker(k).Trajectory.raw)
                                    if isempty(find(isnan(Trial.Marker(k).Trajectory.fill(Trial.Marker(i).Gap(1).frames,1))))
                                        kmarker = [kmarker k];
                                    end
                                end
                            end
                        end
                    end
                end
                
                % If at least 3 markers of the same segment are
                % available, reconstruct the missing marker
                if size(kmarker,2) >= 3
                    X = [];
                    for k = 1:size(kmarker,2)
                        X = [X; Static.Marker(kmarker(k)).Trajectory.fill];
                    end
                    for t = Trial.Marker(i).Gap(1).frames
                        Y = [];
                        for k = 1:size(kmarker,2)
                            Y = [Y; Trial.Marker(kmarker(k)).Trajectory.fill(t,:)];
                        end
                        [R,d,rms] = soder(X,Y);
                        Trial.Marker(i).Trajectory.fill(t,:) = ...
                            Static.Marker(i).Trajectory.fill*R'+d';
                        clear R d;
                    end
                    Trial.Marker(i).Gap(1).reconstruction = 'rigid';
                end
                clear segment;
                
            % Markers related to a curve (semi-landmarks)
            elseif strcmp(Trial.Marker(i).type,'semi-landmark') || ...
                   strcmp(Trial.Marker(i).type,'hybrid-landmark')
                
                % Identify all available markers of the same curve
                % without gap during all frames of the processed gap
                ncurve = Trial.Marker(i).Body.Curve.label;
                kmarker = [];
                if strcmp(ncurve,'none') == 0 % Only for available ncurve
                    for k = 1:size(Trial.Marker,2)
                        if k ~= i
                            if strcmp(Trial.Marker(k).Body.Curve.label,ncurve) == 1
                                if ~isempty(Trial.Marker(k).Trajectory.raw)
                                    if isempty(find(isnan(Trial.Marker(k).Trajectory.fill(Trial.Marker(i).Gap(1).frames,1))))
                                        kmarker = [kmarker k];
                                    end
                                end
                            end
                        end
                    end
                end
                
                % If at least 4 markers of the same curve are available,
                % reconstruct the missing marker
                if size(kmarker,2) >= 4
                    
                    % Set Y level based on the mean of the previous and
                    % next markers
                    prev = 0;
                    next = 0;
                    for k = 1:size(kmarker,2) % distance1 previous and/or next markers
                        if Trial.Marker(kmarker(k)).Body.Curve.index == ...
                                Trial.Marker(i).Body.Curve.index-1
                            prev = k;
                        elseif Trial.Marker(kmarker(k)).Body.Curve.index == ...
                                Trial.Marker(i).Body.Curve.index+1
                            next = k;
                        end
                    end
                    if prev == 0 || next == 0 % distance1 previous and/or next marker not available
                        prev = 0;
                        next = 0;
                        for k = 1:size(kmarker,2) % distance2 previous and/or next markers
                            if Trial.Marker(kmarker(k)).Body.Curve.index == ...
                                    Trial.Marker(i).Body.Curve.index-2
                                prev = k;
                            elseif Trial.Marker(kmarker(k)).Body.Curve.index == ...
                                    Trial.Marker(i).Body.Curve.index+2
                                next = k;
                            end
                        end
                    end
                    if prev ~= 0 && next ~= 0
                        for t = Trial.Marker(i).Gap(1).frames
                            
                            % Store and sort the position of other
                            % semi-landmarks of the related curve
                            temp = [];
                            X = [];
                            for k = 1:size(kmarker,2)
                                temp = [temp; Trial.Marker(kmarker(k)).Body.Curve.index];
                                X    = [X; Trial.Marker(kmarker(k)).Trajectory.fill(t,:)];
                            end
                            [~,I] = sort(temp,'descend');
                            X = X(I,:);
                            clear I;
                            
                            % Define the axis of maximal variation
                            temp = [abs(X(1,1)-X(end,1)) abs(X(1,2)-X(end,2)) abs(X(1,3)-X(end,3))];
                            mVar = find(temp==max(temp));
                            if mVar == 1
                                X1 = X(:,1);
                                % Polynom is plane 1
                                P1 = makima(X1,X(:,2));
                                Y1 = ppval(P1,X1);
                                % Polynom is plane 2
                                P2 = makima(X1,X(:,3));
                                Z1 = ppval(P2,X1);
                                % Reconstruct the missing marker
                                X2 = (Trial.Marker(kmarker(prev)).Trajectory.fill(t,1) + ...
                                    Trial.Marker(kmarker(next)).Trajectory.fill(t,1))/2;
                                Y2 = ppval(P1,X2);
                                Z2 = ppval(P2,X2);
                            elseif mVar == 2
                                Y1 = X(:,2);
                                % Polynom is plane 1
                                P1 = makima(Y1,X(:,1));
                                X1 = ppval(P1,Y1);
                                % Polynom is plane 2
                                P2 = makima(Y1,X(:,3));
                                Z1 = ppval(P2,Y1);
                                % Reconstruct the missing marker
                                Y2 = (Trial.Marker(kmarker(prev)).Trajectory.fill(t,2) + ...
                                    Trial.Marker(kmarker(next)).Trajectory.fill(t,2))/2;
                                X2 = ppval(P1,Y2);
                                Z2 = ppval(P2,Y2);
                            elseif mVar == 3
                                Z1 = X(:,3);
                                % Polynom is plane 1
                                P1 = makima(Z1,X(:,1));
                                X1 = ppval(P1,Z1);
                                % Polynom is plane 2
                                P2 = makima(Z1,X(:,2));
                                Y1 = ppval(P2,Z1);
                                % Reconstruct the missing marker
                                Z2 = (Trial.Marker(kmarker(prev)).Trajectory.fill(t,3) + ...
                                    Trial.Marker(kmarker(next)).Trajectory.fill(t,3))/2;
                                X2 = ppval(P1,Z2);
                                Y2 = ppval(P2,Z2);
                            end
                            Trial.Marker(i).Trajectory.fill(t,:) = ...
                                [X2 Y2 Z2];
                            clear X P1 P2 X2 Y2 Z2;
                        end
                        Trial.Marker(i).Gap(1).reconstruction = 'curve';
                    end
                end
                clear segment;
            end
        end
    end
end

% -------------------------------------------------------------------------
% SMOOTH ALL RESULTING TRAJECTORIES
% -------------------------------------------------------------------------
for i = 1:size(Trial.Marker,2)
    if ~isempty(Trial.Marker(i).Trajectory.fill)
        if ~isempty(Static)
            % Low pass filter (Butterworth 2nd order, [smethod.parameter] Hz)
            if strcmp(smethod.type,'none')    
                Trial.Marker(i).Trajectory.smooth = Trial.Marker(i).Trajectory.fill;
                Trial.Marker(i).smoothing = 'none';
            % Low pass filter (Butterworth 2nd order, [smethod.parameter] Hz)
            elseif strcmp(smethod.type,'butterLow2')                
                [B,A]                             = butter(1,smethod.parameter/(Trial.fmarker/2),'low'); 
                Trial.Marker(i).Trajectory.smooth = filtfilt(B,A,Trial.Marker(i).Trajectory.fill);
                Trial.Marker(i).smoothing = 'butterLow2';
            % Moving average (window of [smethod.parameter] frames)
            elseif strcmp(smethod.type,'movmedian')
                Trial.Marker(i).Trajectory.smooth = smoothdata(Trial.Marker(i).Trajectory.fill,'movmedian',smethod.parameter);
                Trial.Marker(i).smoothing = 'movmedian';
            % Moving average (window of [smethod.parameter] frames)
            elseif strcmp(smethod.type,'movmean')
                Trial.Marker(i).Trajectory.smooth = smoothdata(Trial.Marker(i).Trajectory.fill,'movmean',smethod.parameter);
                Trial.Marker(i).smoothing = 'movmean';
            % Gaussian-weighted moving average (window of [smethod.parameter] frames)
            elseif strcmp(smethod.type,'gaussian')
                Trial.Marker(i).Trajectory.smooth = smoothdata(Trial.Marker(i).Trajectory.fill,'gaussian',smethod.parameter);
                Trial.Marker(i).smoothing = 'gaussian';
            % Robust quadratic regression (window of [smethod.parameter] frames)
            elseif strcmp(smethod.type,'rloess')
                Trial.Marker(i).Trajectory.smooth = smoothdata(Trial.Marker(i).Trajectory.fill,'rloess',smethod.parameter);
                Trial.Marker(i).smoothing = 'rloess';
            % Savitzky-Golay filter (window of [smethod.parameter] frames)
            elseif strcmp(smethod.type,'sgolay')
                Trial.Marker(i).Trajectory.smooth = smoothdata(Trial.Marker(i).Trajectory.fill,'sgolay',smethod.parameter);
                Trial.Marker(i).smoothing = 'sgolay';
            end
        else
            Trial.Marker(i).Trajectory.smooth = Trial.Marker(i).Trajectory.fill;
        end
        % Modify the ICS (Y vertical)
        Trial.Marker(i).Trajectory.smooth = [Trial.Marker(i).Trajectory.smooth(:,1) ...
                                             Trial.Marker(i).Trajectory.smooth(:,3) ...
                                             -Trial.Marker(i).Trajectory.smooth(:,2)];
    else
        Trial.Marker(i).Trajectory.smooth = [];
    end
end

% for i = 1:size(Trial.Marker,2)
%     if ~isempty(Trial.Marker(i).trajectory)
%         if ~isempty(Trial.Marker(i).Gap)
%             for j = 1:size(Trial.Marker(i).Gap,2)
%                 if strcmp(Trial.Marker(i).Gap(j).type,'Short gap')
%                     temp = PredictMissingMarkers([fullMarker Trial.Marker(i).trajectory]);
%                     Trial.Marker(i).trajectory(Trial.Marker(i).Gap(j).frames,:) = temp(Trial.Marker(i).Gap(j).frames,end-2:end);
%                     Trial.Marker(i).Gap(j).reconstruction = 'PCA';
%                     clear temp;
%                 end
%             end
%         end
%     end
% end

% Step 2: Reconstruct long gaps >= gap threshold with a rigid body or curve 
%         approach (also applied on missing markers)
% if ~isempty(Static)
%     for i = 1:size(Trial.Marker,2)
%         if ~isempty(Trial.Marker(i).Gap)
%             for j = 1:size(Trial.Marker(i).Gap,2)
%                 if strcmp(Trial.Marker(i).Gap(j).type,'Long gap') == 1 || ...
%                    strcmp(Trial.Marker(i).Gap(j).type,'Missing trajectory') == 1
% 
%                     % Markers related to a rigid body (landmarks and hybrid-landmarks)
%                     if strcmp(Trial.Marker(i).type,'landmark') || ...
%                        strcmp(Trial.Marker(i).type,'hybrid-landmark') || ...
%                        strcmp(Trial.Marker(i).type,'technical')
% 
%                         % Identify all available markers of the same segment 
%                         % without gap during all frames of the processed gap
%                         nsegment = Trial.Marker(i).Body.Segment.label;
%                         kmarker = [];
%                         if strcmp(nsegment,'none') == 0 % Only for available segments                            
%                             for k = 1:size(Trial.Marker,2)
%                                 if k ~= i
%                                     if strcmp(Trial.Marker(k).Body.Segment.label,nsegment) == 1
%                                         if ~isempty(Trial.Marker(k).trajectory)
%                                             if isempty(find(isnan(Trial.Marker(k).trajectory(Trial.Marker(i).Gap(j).frames,1))))
%                                                 kmarker = [kmarker k];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
% 
%                         % If at least 3 markers of the same segment are 
%                         % available, reconstruct the missing marker
%                         if size(kmarker,2) >= 3
%                             X = [];
%                             for k = 1:size(kmarker,2)
%                                 X = [X; Static.Marker(kmarker(k)).trajectory];
%                             end
%                             for t = Trial.Marker(i).Gap(j).frames
%                                 Y = [];
%                                 for k = 1:size(kmarker,2)
%                                     Y = [Y; Trial.Marker(kmarker(k)).trajectory(t,:)];
%                                 end
%                                 [R,d,rms] = soder(X,Y); 
%                                 Trial.Marker(i).trajectory(t,:) = ...
%                                     Static.Marker(i).trajectory*R'+d';
%                                 clear R d;
%                             end                            
%                             Trial.Marker(i).Gap(j).reconstruction = 'Rigid body';
%                             disp(['      > ',Trial.Marker(i).label,' gap',num2str(j),' filled by rigid body']);
%                         else
%                             disp(['      > ',Trial.Marker(i).label,' gap',num2str(j),' not filled (not enough markers on the rigid segment)']);
%                         end
%                         clear segment;
%                     
%                     % Markers related to a curve (semi-landmarks)
%                     elseif strcmp(Trial.Marker(i).type,'semi-landmark') || ...
%                            strcmp(Trial.Marker(i).type,'hybrid-landmark')
% 
%                         % Identify all available markers of the same curve 
%                         % without gap during all frames of the processed gap
%                         ncurve = Trial.Marker(i).Body.Curve.label;
%                         kmarker = [];
%                         if strcmp(ncurve,'none') == 0 % Only for available ncurve                            
%                             for k = 1:size(Trial.Marker,2)
%                                 if k ~= i
%                                     if strcmp(Trial.Marker(k).Body.Curve.label,ncurve) == 1
%                                         if ~isempty(Trial.Marker(k).trajectory)
%                                             if isempty(find(isnan(Trial.Marker(k).trajectory(Trial.Marker(i).Gap(j).frames,1))))
%                                                 kmarker = [kmarker k];
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     
%                         % If at least 4 markers of the same curve are available,
%                         % reconstruct the missing marker
%                         if size(kmarker,2) >= 4
% 
%                             % Set Y level based on the mean of the previous and
%                             % next markers
%                             prev = 0;
%                             next = 0;
%                             for k = 1:size(kmarker,2) % distance1 previous and/or next markers
%                                 if Trial.Marker(kmarker(k)).Body.Curve.index == ...
%                                    Trial.Marker(i).Body.Curve.index-1
%                                     prev = k;
%                                 elseif Trial.Marker(kmarker(k)).Body.Curve.index == ...
%                                        Trial.Marker(i).Body.Curve.index+1
%                                     next = k;
%                                 end
%                             end
%                             if prev == 0 || next == 0 % distance1 previous and/or next marker not available
%                                 prev = 0;
%                                 next = 0;
%                                 for k = 1:size(kmarker,2) % distance2 previous and/or next markers
%                                     if Trial.Marker(kmarker(k)).Body.Curve.index == ...
%                                        Trial.Marker(i).Body.Curve.index-2
%                                         prev = k;
%                                     elseif Trial.Marker(kmarker(k)).Body.Curve.index == ...
%                                            Trial.Marker(i).Body.Curve.index+2
%                                         next = k;
%                                     end
%                                 end
%                             end
%                             if prev == 0 || next == 0
%                                 disp(['      > ',Trial.Marker(i).label,' large gaps not filled (not enough close markers on the related curve)']);
%                             else              
%                                 for t = Trial.Marker(i).Gap(j).frames
%                                     % Store and sort the position of other
%                                     % semi-landmarks of the related curve
%                                     temp = [];
%                                     X = [];
%                                     for k = 1:size(kmarker,2)
%                                         temp = [temp; Trial.Marker(kmarker(k)).Body.Curve.index];
%                                         X    = [X; Trial.Marker(kmarker(k)).trajectory(t,:)];
%                                     end
%                                     [~,I] = sort(temp,'descend');
%                                     X = X(I,:);
%                                     clear I;
%                                     % Define the axis of maximal variation
%                                     temp = [abs(X(1,1)-X(end,1)) abs(X(1,2)-X(end,2)) abs(X(1,3)-X(end,3))];
%                                     mVar = find(temp==max(temp));
%                                     if mVar == 1
%                                         X1 = X(:,1);
%                                         % Polynom is plane 1
%                                         P1 = makima(X1,X(:,2));
%                                         Y1 = ppval(P1,X1);
%                                         % Polynom is plane 2
%                                         P2 = makima(X1,X(:,3));
%                                         Z1 = ppval(P2,X1);
%                                         % Reconstruct the missing marker
%                                         X2 = (Trial.Marker(kmarker(prev)).trajectory(t,1) + ...
%                                               Trial.Marker(kmarker(next)).trajectory(t,1))/2;
%                                         Y2 = ppval(P1,X2);
%                                         Z2 = ppval(P2,X2);
%                                     elseif mVar == 2
%                                         Y1 = X(:,2);
%                                         % Polynom is plane 1
%                                         P1 = makima(Y1,X(:,1));
%                                         X1 = ppval(P1,Y1);
%                                         % Polynom is plane 2
%                                         P2 = makima(Y1,X(:,3));
%                                         Z1 = ppval(P2,Y1);
%                                         % Reconstruct the missing marker
%                                         Y2 = (Trial.Marker(kmarker(prev)).trajectory(t,2) + ...
%                                               Trial.Marker(kmarker(next)).trajectory(t,2))/2;
%                                         X2 = ppval(P1,Y2);
%                                         Z2 = ppval(P2,Y2);
%                                     elseif mVar == 3
%                                         Z1 = X(:,3);
%                                         % Polynom is plane 1
%                                         P1 = makima(Z1,X(:,1));
%                                         X1 = ppval(P1,Z1);
%                                         % Polynom is plane 2
%                                         P2 = makima(Z1,X(:,2));
%                                         Y1 = ppval(P2,Z1);
%                                         % Reconstruct the missing marker
%                                         Z2 = (Trial.Marker(kmarker(prev)).trajectory(t,3) + ...
%                                               Trial.Marker(kmarker(next)).trajectory(t,3))/2;
%                                         X2 = ppval(P1,Z2);
%                                         Y2 = ppval(P2,Z2);
%                                     end
%                                     Trial.Marker(i).trajectory(t,:) = ...
%                                         [X2 Y2 Z2];
%                                     clear X P1 P2 X2 Y2 Z2;
%                                 end                                
%                                 Trial.Marker(i).Gap(j).reconstruction = 'Curve';
%                                 disp(['      > ',Trial.Marker(i).label,' large gaps filled (along a curve)']);
%                             end
%                         else
%                             disp(['      > ',Trial.Marker(i).label,' large gaps not filled (not enough markers on the related curve)']);
%                         end
%                         clear segment;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % Step 3: Smooth marker trajectories (or keep only the mean position for 
% %         static)
% for i = 1:size(Trial.Marker,2)
%     if ~isempty(Trial.Marker(i).Trajectory.raw)
%         if ~isempty(Static)
%             % Low pass filter (Butterworth 2nd order, 6 Hz)
%             [B,A]                      = butter(1,6/(Trial.fmarker/2),'low'); 
%             Trial.Marker(i).trajectory = filtfilt(B,A,Trial.Marker(i).trajectory);
%         elseif isempty(Static)
%             Trial.Marker(i).trajectory = mean(Trial.Marker(i).trajectory,1);
%             Trial.n0                   = 1;
%             Trial.n1                   = 1;
%         end
%         % Modify the ICS (Y vertical)
%         if ~isempty(Trial.Marker(i).trajectory)
%             Trial.Marker(i).trajectory = [Trial.Marker(i).trajectory(:,1) ...
%                                           Trial.Marker(i).trajectory(:,3) ...
%                                           -Trial.Marker(i).trajectory(:,2)];
%         end
%     end
% end

% 
% fmethod.gapThreshold = 100*Trial.fmarker/100; % frames (1 second allowed)
% fullMarkers  = []; % Merge of markers without any gap (needed for PCA reconstruction)
% nMarker      = fieldnames(Trial.Marker);
% for i = 1:size(nMarker,1)
%     if isempty(Trial.Marker.(nMarker{i}).trajectory)
%         disp(['      > ',nMarker{i},' is missing']);
%     elseif ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%         % Initialise gap field
%         igap                           = [];
%         Trial.Marker.(nMarker{i}).gap  = {}; % Column 1: 
%         Trial.Marker.(nMarker{i}).igap = [];
%         % Replace [0 0 0] by NaN
%         for j = 1:Trial.n1
%             if Trial.Marker.(nMarker{i}).trajectory(j,:) == [0 0 0]
%                Trial.Marker.(nMarker{i}).trajectory(j,:) = nan(1,3);
%             end
%         end
%         % Find gaps
%         gap   = [];
%         start = 0;
%         stop  = 0;
%         k     = 0;
%         for j = 1:Trial.n1-1
%             if isnan(Trial.Marker.(nMarker{i}).trajectory(j,:))
%                 if start == 0
%                     start = j;
%                 end
%                 if ~isnan(Trial.Marker.(nMarker{i}).trajectory(j+1,:))
%                     if start ~= 0
%                         stop   = j;
%                         k      = k+1;    
%                         gap{k} = start:stop;                    
%                         start  = 0;
%                         stop   = 0;
%                         % Group gaps < or > the gap threshold
%                         if size(gap{k},2) >= fmethod.gapThreshold
%                             igap(k) = 1;
%                         else
%                             igap(k) = 2;
%                         end
%                     end
%                 elseif j+1 == Trial.n1
%                     if isnan(Trial.Marker.(nMarker{i}).trajectory(j+1,:))
%                         if start ~= 0
%                             stop   = j+1;
%                             k      = k+1;    
%                             gap{k} = start:stop;                    
%                             start  = 0;
%                             stop   = 0;
%                             % Group gaps < or > the gap threshold
%                             if size(gap{k},2) >= fmethod.gapThreshold
%                                 igap(k) = 1;
%                             else
%                                 igap(k) = 2;
%                             end
%                         end   
%                     end
%                 end
%             end
%         end
%         if k == 0 % No gap found
%             fullMarkers = [fullMarkers Trial.Marker.(nMarker{i}).trajectory];
%         end
%         Trial.Marker.(nMarker{i}).gap = gap;
%         Trial.Marker.(nMarker{i}).igap = igap;
%     end
% end
% 
% % -------------------------------------------------------------------------
% % PROCESS MARKER TRAJECTORIES
% % -------------------------------------------------------------------------
% 
% % Step 1: Reconstruct only missing trajectories of markers having only gaps
% %         < gap threshold with a PCA fmethod.type based on full marker
% %         trajectories (no gap)
% for i = 1:size(nMarker,1)
%     if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%         if ~isempty(Trial.Marker.(nMarker{i}).igap)
%             if isempty(find(Trial.Marker.(nMarker{i}).igap == 1))
%                 temp = PredictMissingMarkers([fullMarkers Trial.Marker.(nMarker{i}).trajectory]);
%                 Trial.Marker.(nMarker{i}).trajectory = temp(:,end-2:end);
%                 clear temp;
%             end
%         end
%     end
% end
% 
% % Step 2: Reconstruct only missing trajectories of markers having only gaps
% %         >= gap threshold with a rigid body or curve approach (with full
% %         markers only)
% for i = 1:size(nMarker,1)
%     if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%         if ~isempty(Trial.Marker.(nMarker{i}).igap)
%             if ~isempty(find(Trial.Marker.(nMarker{i}).igap == 1))
% 
%                 % Markers related to a rigid body (landmarks and hybrid-landmarks)
%                 if strcmp(Trial.Marker.(nMarker{i}).type,'landmark') || ...
%                    strcmp(Trial.Marker.(nMarker{i}).type,'hybrid-landmark')
%                     
%                     % Define related segment
%                     nsegment = Trial.Marker.(nMarker{i}).segment;
%                
%                     for k = 1:size(Trial.Marker.(nMarker{i}).gap,2)
% 
%                         % Identify available full markers of the same segment
%                         jmarker = [];
%                         if strcmp(nsegment,'none') == 0 % Only for available segments
%                             for j = 1:size(nMarker,1)
%                                 if strcmp(Trial.Marker.(nMarker{j}).segment,nsegment) == 1 % Marker of the same segment
%                                     if isfield(Trial.Marker.(nMarker{j}),'igap')
%                                         if isempty(find(Trial.Marker.(nMarker{j}).igap == 1))
%                                             jmarker = [jmarker j];
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                         
%                         % If at least 3 markers of the same segment are 
%                         % available, reconstruct the missing marker
%                         if size(jmarker,2) >= 3
%                             X = [];
%                             for j = 1:size(jmarker,2)
%                                 X = [X; Static.Marker.(nMarker{jmarker(j)}).trajectory];
%                             end
%                             for t = Trial.Marker.(nMarker{i}).gap{k}(1):Trial.Marker.(nMarker{i}).gap{k}(end)
%                                 Y = [];
%                                 for j = 1:size(jmarker,2)
%                                     Y = [Y; Trial.Marker.(nMarker{jmarker(j)}).trajectory(t,:)];
%                                 end
%                                 [R,d,rms] = soder(X,Y); 
%                                 Trial.Marker.(nMarker{i}).trajectory(t,:) = ...
%                                     Static.Marker.(nMarker{i}).trajectory*R'+d';
%                                 clear Y R d;
%                             end
%                             disp(['      > ',nMarker{i},' gap',num2str(k),' filled by rigid body']);
%                         else
%                             disp(['      > ',nMarker{i},' gap',num2str(k),' not filled (not enough markers on the rigid segment)']);
%                         end
%                         clear segment;
%                     end
% 
%                 % Markers related to a curve (semi-landmarks)
%                 elseif strcmp(Trial.Marker.(nMarker{i}).type,'semi-landmark') || ...
%                        strcmp(Trial.Marker.(nMarker{i}).type,'hybrid-landmark')
% 
% %                     % Define related curve
% %                     curve = Trial.Marker.(nMarker{i}).curve{1};
% %                     % Identify available markers of the same curve
% %                     jmarker = [];
% %                     if strcmp(curve,'none') == 0 % Only for available curves
% %                         for j = 1:size(nMarker,1)
% %                             if j ~= i % Other marker than the one missing
% %                                 if strcmp(Trial.Marker.(nMarker{j}).curve{1},curve) == 1 % Marker of the same curve
% %                                     if ~isempty(Trial.Marker.(nMarker{j}).trajectory)
% %                                         if ~isempty(Trial.Marker.(nMarker{j}).trajectory) % Non-empty marker
% %                                             jmarker = [jmarker j];
% %                                         end
% %                                     end
% %                                 end
% %                             end
% %                         end
% %                     end            
% %                     % If at least 4 markers of the same curve are available,
% %                     % reconstruct the missing marker
% %                     for k = 1:size(Trial.Marker.(nMarker{i}).gap,2)
% %                         if size(jmarker,2) >= 4
% %     
% %                             % Set Y level based on the mean of the previous and
% %                             % next markers
% %                             prev = 0;
% %                             next = 0;
% %                             for j = 1:size(jmarker,2) % distance1 previous and/or next markers
% %                                 if Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
% %                                    Trial.Marker.(nMarker{i}).curve{2}-1
% %                                     prev = j;
% %                                 elseif Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
% %                                        Trial.Marker.(nMarker{i}).curve{2}+1
% %                                     next = j;
% %                                 end
% %                             end
% %                             if prev == 0 || next == 0 % distance1 previous and/or next marker not available
% %                                 prev = 0;
% %                                 next = 0;
% %                                 for j = 1:size(jmarker,2) % distance2 previous and/or next markers
% %                                     if Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
% %                                        Trial.Marker.(nMarker{i}).curve{2}-2
% %                                         prev = j;
% %                                     elseif Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
% %                                            Trial.Marker.(nMarker{i}).curve{2}+2
% %                                         next = j;
% %                                     end
% %                                 end
% %                             end
% %                             if prev == 0 || next == 0
% %                                 disp(['      > ',nMarker{i},' large gaps not filled (not enough close markers on the related curve)']);
% %                             else              
% %                                 for t = Trial.Marker.(nMarker{i}).gap{k}(1):Trial.Marker.(nMarker{i}).gap{k}(end)
% %                                     % Store and sort the position of other
% %                                     % semi-landmarks of the related curve
% %                                     temp = [];
% %                                     X = [];
% %                                     for j = 1:size(jmarker,2)
% %                                         temp = [temp; Trial.Marker.(nMarker{jmarker(j)}).curve{2}];
% %                                         X    = [X; Trial.Marker.(nMarker{jmarker(j)}).trajectory(t,:)];
% %                                     end
% %                                     [~,I] = sort(temp,'descend');
% %                                     X = X(I,:);
% %                                     clear I;
% %                                     % Define the axis of maximal variation
% %                                     temp = [abs(X(1,1)-X(end,1)) abs(X(1,2)-X(end,2)) abs(X(1,3)-X(end,3))];
% %                                     mVar = find(temp==max(temp));
% %                                     if mVar == 1
% %                                         X1 = X(:,1);
% %                                         % Polynom is plane 1
% %                                         P1 = makima(X1,X(:,2));
% %                                         Y1 = ppval(P1,X1);
% %                                         % Polynom is plane 2
% %                                         P2 = makima(X1,X(:,3));
% %                                         Z1 = ppval(P2,X1);
% %                                         % Reconstruct the missing marker
% %                                         X2 = (Trial.Marker.(nMarker{jmarker(prev)}).trajectory(t,1) + ...
% %                                               Trial.Marker.(nMarker{jmarker(next)}).trajectory(t,1))/2;
% %                                         Y2 = ppval(P1,X2);
% %                                         Z2 = ppval(P2,X2);
% %                                     elseif mVar == 2
% %                                         Y1 = X(:,2);
% %                                         % Polynom is plane 1
% %                                         P1 = makima(Y1,X(:,1));
% %                                         X1 = ppval(P1,Y1);
% %                                         % Polynom is plane 2
% %                                         P2 = makima(Y1,X(:,3));
% %                                         Z1 = ppval(P2,Y1);
% %                                         % Reconstruct the missing marker
% %                                         Y2 = (Trial.Marker.(nMarker{jmarker(prev)}).trajectory(t,2) + ...
% %                                               Trial.Marker.(nMarker{jmarker(next)}).trajectory(t,2))/2;
% %                                         X2 = ppval(P1,Y2);
% %                                         Z2 = ppval(P2,Y2);
% %                                     elseif mVar == 3
% %                                         Z1 = X(:,3);
% %                                         % Polynom is plane 1
% %                                         P1 = makima(Z1,X(:,1));
% %                                         X1 = ppval(P1,Z1);
% %                                         % Polynom is plane 2
% %                                         P2 = makima(Z1,X(:,2));
% %                                         Y1 = ppval(P2,Z1);
% %                                         % Reconstruct the missing marker
% %                                         Z2 = (Trial.Marker.(nMarker{jmarker(prev)}).trajectory(t,3) + ...
% %                                               Trial.Marker.(nMarker{jmarker(next)}).trajectory(t,3))/2;
% %                                         X2 = ppval(P1,Z2);
% %                                         Y2 = ppval(P2,Z2);
% %                                     end
% %                                     Trial.Marker.(nMarker{i}).trajectory(t,:) = ...
% %                                         [X2 Y2 Z2];
% %                                     clear X P1 P2 X2 Y2 Z2;
% %                                 end
% %                                 disp(['      > ',nMarker{i},' large gaps filled (along a curve)']);
% %                             end
% %                         else
% %                             disp(['      > ',nMarker{i},' large gaps not filled (not enough markers on the related curve)']);
% %                         end
% %                         clear segment;
% %                     end
% 
%                 end
%             end
%         end
%     end
% end
% 
% % Step 2: Smooth marker trajectories (or keep only the mean position for 
% %         static)
% for i = 1:size(nMarker,1)
%     if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%         if ~isempty(Static)
% %             % Bandstop filter (Butterworth 2nd order, exclude 4-40 Hz)
% %             [B,A]                                = butter(4,[4/(Trial.fmarker/2) 40/(Trial.fmarker/2)],'stop'); 
% %             Trial.Marker.(nMarker{i}).trajectory = filtfilt(B,A,Trial.Marker.(nMarker{i}).trajectory);
%             % Low pass filter (Butterworth 2nd order, 6 Hz)
%             [B,A]                                = butter(1,6/(Trial.fmarker/2),'low'); 
% %             Trial.Marker.(nMarker{i}).trajectory = filtfilt(B,A,Trial.Marker.(nMarker{i}).trajectory);
%         elseif isempty(Static)
%             Trial.Marker.(nMarker{i}).trajectory = mean(Trial.Marker.(nMarker{i}).trajectory,1);
%             Trial.n0                             = 1;
%             Trial.n1                             = 1;
%         end
%         % Convert mm to m
%         Trial.Marker.(nMarker{i}).trajectory = Trial.Marker.(nMarker{i}).trajectory*1e-3;
%         % Modify the ICS (Y vertical)
%         if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%             Trial.Marker.(nMarker{i}).trajectory = [Trial.Marker.(nMarker{i}).trajectory(:,1) ...
%                                                     Trial.Marker.(nMarker{i}).trajectory(:,3) ...
%                                                     -Trial.Marker.(nMarker{i}).trajectory(:,2)];
%         end
%     end
% end

% % -------------------------------------------------------------------------
% % PROCESS AVAILABLE STATIC AND TRIAL MARKER TRAJECTORIES
% % -------------------------------------------------------------------------
% imarker = [];
% nMarker = fieldnames(Trial.Marker);
% for i = 1:size(nMarker,1)
%     if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%         % Replace [0 0 0] by NaN
%         for j = 1:Trial.n1
%             if Trial.Marker.(nMarker{i}).trajectory(j,:) == [0 0 0]
%                Trial.Marker.(nMarker{i}).trajectory(j,:) = nan(1,3);
%             end
%         end    
% %         % Find marker trajectories where too many frames are missing
% %         threshold = 20; % At least [threshold]% of the marker trajectory is missing
% %         if ~isempty(Static)
% %             temp = isnan(Trial.Marker.(nMarker{i}).trajectory);
% %             if ceil(sum(temp(:,1))*100/Trial.n1) >= threshold
% %                 imarker = [imarker i];
% %             end
% %         end
%         % Find marker trajectories where gaps are too long
%         threshold = 300*Trial.fmarker/100; % Number of consecutive frames (initially considered at 100 Hz)
%         gap    = 0;
%         maxGap = 0;
%         if isnan(Trial.Marker.(nMarker{i}).trajectory(1,1)) % First frame
%             gap = gap+1;
%             maxGap = gap;
%         end
%         for j = 2:Trial.n1
%             if isnan(Trial.Marker.(nMarker{i}).trajectory(j-1,1)) && ...
%                isnan(Trial.Marker.(nMarker{i}).trajectory(j,1)) % Tested on dimension 1 only
%                 gap = gap+1;
%                 if gap > maxGap
%                     maxGap = gap;
%                 end
%             else
%                 gap = 0;
%             end
%         end
%         if maxGap > threshold
%             imarker = [imarker i];
%         end
%     end
% end
% % Clear marker trajectories when too many frames are missing
% for i = 1:size(imarker,2)
%     % Clear marker
%     Trial.Marker.(nMarker{imarker(i)}).trajectory = [];
%     disp(['      > ',nMarker{imarker(i)},' cleared (gaps >',num2str(threshold),' frames)']);
% end
% clear imarker;
% % Merge marker data for interpolation 
% tMarker = [];
% for i = 1:size(nMarker,1) 
%     if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%         tMarker = [tMarker Trial.Marker.(nMarker{i}).trajectory];
%     end
% end
% % Interpolate markers trajectory
% tMarker = PredictMissingMarkers(tMarker);
% j       = 1;
% % Smooth marker trajectories (or keep only the mean position for static)
% for i = 1:size(nMarker,1)
%     if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%         if ~isempty(Static)
%             % Extract interpolated data        
%             Trial.Marker.(nMarker{i}).trajectory = tMarker(:,(3*j)-2:3*j);
%             j = j+1;
% %             % Bandstop filter (Butterworth 2nd order, exclude 4-40 Hz)
% %             [B,A]                                = butter(4,[4/(Trial.fmarker/2) 40/(Trial.fmarker/2)],'stop'); 
% %             Trial.Marker.(nMarker{i}).trajectory = filtfilt(B,A,Trial.Marker.(nMarker{i}).trajectory);
%             % Low pass filter (Butterworth 2nd order, 6 Hz)
%             [B,A]                                = butter(1,6/(Trial.fmarker/2),'low'); 
%             Trial.Marker.(nMarker{i}).trajectory = filtfilt(B,A,Trial.Marker.(nMarker{i}).trajectory);
%         elseif isempty(Static)
%             Trial.Marker.(nMarker{i}).trajectory = mean(Trial.Marker.(nMarker{i}).trajectory,1);
%             Trial.n0                             = 1;
%             Trial.n1                             = 1;
%         end
%         % Convert mm to m
%         Trial.Marker.(nMarker{i}).trajectory = Trial.Marker.(nMarker{i}).trajectory*1e-3;
%         % Modify the ICS (Y vertical)
%         if ~isempty(Trial.Marker.(nMarker{i}).trajectory)
%             Trial.Marker.(nMarker{i}).trajectory = [Trial.Marker.(nMarker{i}).trajectory(:,1) ...
%                                                     Trial.Marker.(nMarker{i}).trajectory(:,3) ...
%                                                     -Trial.Marker.(nMarker{i}).trajectory(:,2)];
%         end
%     end
% end
% 
% % -------------------------------------------------------------------------
% % RECONSTRUCT UNAVAILABLE TRIAL MARKER TRAJECTORIES
% % -------------------------------------------------------------------------
% if ~isempty(Static)
%     
%     % Find missing markers
%     imarker = [];
%     for i = 1:size(nMarker,1)
%         if isempty(Trial.Marker.(nMarker{i}).trajectory)
%             imarker = [imarker i];
%         end
%     end
%     
%     for i = 1:size(imarker,2)
%         
%         % Markers related to a rigid body (landmarks and hybrid-landmarks)
%         if strcmp(Trial.Marker.(nMarker{imarker(i)}).type,'landmark') || ...
%            strcmp(Trial.Marker.(nMarker{imarker(i)}).type,'hybrid-landmark')
%         
%             % Define related segment
%             segment = Trial.Marker.(nMarker{imarker(i)}).segment;
%             % Identify available markers of the same segment
%             jmarker = [];
%             if strcmp(segment,'none') == 0 % Only for available segments
%                 for j = 1:size(nMarker,1)
%                     if j ~= i % Other marker than the one missing
%                         if strcmp(Trial.Marker.(nMarker{j}).segment,segment) == 1 % Marker of the same segment
%                             if ~isempty(Trial.Marker.(nMarker{j}).trajectory) % Non-empty marker
%                                 jmarker = [jmarker j];
%                             end
%                         end
%                     end
%                 end
%             end
%             % If at least 3 markers of the same segment are available,
%             % reconstruct the missing marker
%             if size(jmarker,2) >= 3
%                 X = [];
%                 for j = 1:size(jmarker,2)
%                     X = [X; Static.Marker.(nMarker{jmarker(j)}).trajectory];
%                 end
%                 for t = 1:Trial.n1
%                     Y = [];
%                     for j = 1:size(jmarker,2)
%                         Y = [Y; Trial.Marker.(nMarker{jmarker(j)}).trajectory(t,:)];
%                     end
%                     [R,d,rms] = soder(X,Y); 
%                     Trial.Marker.(nMarker{imarker(i)}).trajectory(t,:) = ...
%                         Static.Marker.(nMarker{imarker(i)}).trajectory*R'+d';
%                     clear Y R d;
%                 end
%                 disp(['      > ',nMarker{imarker(i)},' reconstructed (rigid body)']);
%                 clear X;
%             else
%                 disp(['      > ',nMarker{imarker(i)},' not reconstructed (not enough markers on the rigid segment)']);
%             end
%             clear segment;
%             
%         % Markers related to a curve (semi-landmarks)
%         elseif strcmp(Trial.Marker.(nMarker{imarker(i)}).type,'semi-landmark') || ...
%                strcmp(Trial.Marker.(nMarker{imarker(i)}).type,'hybrid-landmark')
%             
%             % Define related curve
%             curve = Trial.Marker.(nMarker{imarker(i)}).curve{1};
%             % Identify available markers of the same curve
%             jmarker = [];
%             if strcmp(curve,'none') == 0 % Only for available curves
%                 for j = 1:size(nMarker,1)
%                     if j ~= i % Other marker than the one missing
%                         if strcmp(Trial.Marker.(nMarker{j}).curve{1},curve) == 1 % Marker of the same curve
%                             if ~isempty(Trial.Marker.(nMarker{j}).trajectory) % Non-empty marker
%                                 jmarker = [jmarker j];
%                             end
%                         end
%                     end
%                 end
%             end            
%             % If at least 4 markers of the same curve are available,
%             % reconstruct the missing marker
%             if size(jmarker,2) >= 4
% 
%                 % Set Y level based on the mean of the previous and
%                 % next markers
%                 prev = 0;
%                 next = 0;
%                 for j = 1:size(jmarker,2) % distance1 previous and/or next markers
%                     if Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
%                        Trial.Marker.(nMarker{imarker(i)}).curve{2}-1
%                         prev = j;
%                     elseif Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
%                            Trial.Marker.(nMarker{imarker(i)}).curve{2}+1
%                         next = j;
%                     end
%                 end
%                 if prev == 0 || next == 0 % distance1 previous and/or next marker not available
%                     prev = 0;
%                     next = 0;
%                     for j = 1:size(jmarker,2) % distance2 previous and/or next markers
%                         if Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
%                            Trial.Marker.(nMarker{imarker(i)}).curve{2}-2
%                             prev = j;
%                         elseif Trial.Marker.(nMarker{jmarker(j)}).curve{2} == ...
%                                Trial.Marker.(nMarker{imarker(i)}).curve{2}+2
%                             next = j;
%                         end
%                     end
%                 end
%                 if prev == 0 || next == 0
%                     disp(['      > ',nMarker{imarker(i)},' not reconstructed (not enough close markers on the related curve)']);
%                 else              
%                     for t = 1:Trial.n1
%                         % Store and sort the position of other
%                         % semi-landmarks of the related curve
%                         temp = [];
%                         X = [];
%                         for j = 1:size(jmarker,2)
%                             temp = [temp; Trial.Marker.(nMarker{jmarker(j)}).curve{2}];
%                             X    = [X; Trial.Marker.(nMarker{jmarker(j)}).trajectory(t,:)];
%                         end
%                         [~,I] = sort(temp,'descend');
%                         X = X(I,:);
%                         clear I;
%                         % Define the axis of maximal variation
%                         temp = [abs(X(1,1)-X(end,1)) abs(X(1,2)-X(end,2)) abs(X(1,3)-X(end,3))];
%                         mVar = find(temp==max(temp));
%                         if mVar == 1
%                             X1 = X(:,1);
%                             % Polynom is plane 1
%                             P1 = makima(X1,X(:,2));
%                             Y1 = ppval(P1,X1);
%                             % Polynom is plane 2
%                             P2 = makima(X1,X(:,3));
%                             Z1 = ppval(P2,X1);
%                             % Reconstruct the missing marker
%                             X2 = (Trial.Marker.(nMarker{jmarker(prev)}).trajectory(t,1) + ...
%                                   Trial.Marker.(nMarker{jmarker(next)}).trajectory(t,1))/2;
%                             Y2 = ppval(P1,X2);
%                             Z2 = ppval(P2,X2);
%                         elseif mVar == 2
%                             Y1 = X(:,2);
%                             % Polynom is plane 1
%                             P1 = makima(Y1,X(:,1));
%                             X1 = ppval(P1,Y1);
%                             % Polynom is plane 2
%                             P2 = makima(Y1,X(:,3));
%                             Z1 = ppval(P2,Y1);
%                             % Reconstruct the missing marker
%                             Y2 = (Trial.Marker.(nMarker{jmarker(prev)}).trajectory(t,2) + ...
%                                   Trial.Marker.(nMarker{jmarker(next)}).trajectory(t,2))/2;
%                             X2 = ppval(P1,Y2);
%                             Z2 = ppval(P2,Y2);
%                         elseif mVar == 3
%                             Z1 = X(:,3);
%                             % Polynom is plane 1
%                             P1 = makima(Z1,X(:,1));
%                             X1 = ppval(P1,Z1);
%                             % Polynom is plane 2
%                             P2 = makima(Z1,X(:,2));
%                             Y1 = ppval(P2,Z1);
%                             % Reconstruct the missing marker
%                             Z2 = (Trial.Marker.(nMarker{jmarker(prev)}).trajectory(t,3) + ...
%                                   Trial.Marker.(nMarker{jmarker(next)}).trajectory(t,3))/2;
%                             X2 = ppval(P1,Z2);
%                             Y2 = ppval(P2,Z2);
%                         end
%                         Trial.Marker.(nMarker{imarker(i)}).trajectory(t,:) = ...
%                             [X2 Y2 Z2];
%                         clear X P1 P2 X2 Y2 Z2;
%                     end
%                     disp(['      > ',nMarker{imarker(i)},' reconstructed (along a curve)']);
%                 end
%             else
%                 disp(['      > ',nMarker{imarker(i)},' not reconstructed (not enough markers on the related curve)']);
%             end
%             clear segment;
%             
%         end            
%     end    
% end