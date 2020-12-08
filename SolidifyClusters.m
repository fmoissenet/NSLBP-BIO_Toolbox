% Author       : F. Moissenet
%                Kinesiology Laboratory (K-LAB)
%                University of Geneva
%                https://www.unige.ch/medecine/kinesiology
% License      : Creative Commons Attribution-NonCommercial 4.0 International License 
%                https://creativecommons.org/licenses/by-nc/4.0/legalcode
% Source code  : https://github.com/fmoissenet/Roboshoulder_Toolbox
% Reference    : To be defined
% Date         : March 2020
% -------------------------------------------------------------------------
% Description  : This routine aims to ...
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

function [M1,M2,M3,M4] = SolidifyClusters(M1,M2,M3,M4)

% -------------------------------------------------------------------------
% DETERMINATION OF A SOLID TRIANGLE
% -------------------------------------------------------------------------

% Initialise variables
marker = NaN(size(M1,1),3,4);
vector = NaN(size(M1,1),3,6);
theta  = NaN(size(M1,1),1,9);

% Store cluster marker trajectories
marker(:,:,1) = M1;
marker(:,:,2) = M2;
marker(:,:,3) = M3;
marker(:,:,4) = M4;

% Set all related vectors
vector(:,:,1) = marker(:,:,2)-marker(:,:,1);
vector(:,:,2) = marker(:,:,3)-marker(:,:,1);
vector(:,:,3) = marker(:,:,3)-marker(:,:,2);
vector(:,:,4) = marker(:,:,4)-marker(:,:,1);
vector(:,:,5) = marker(:,:,4)-marker(:,:,3);
vector(:,:,6) = marker(:,:,4)-marker(:,:,2);

% Set all related angles
for i = 1:size(M1,1)
    theta(i,:,1)  = atan2(norm(cross(vector(i,:,1),vector(i,:,2))),dot(vector(i,:,1),vector(i,:,2)));
    theta(i,:,2)  = atan2(norm(cross(vector(i,:,1),vector(i,:,3))),dot(vector(i,:,1),vector(i,:,3)));
    theta(i,:,3)  = atan2(norm(cross(vector(i,:,2),vector(i,:,3))),dot(vector(i,:,2),vector(i,:,3)));
    theta(i,:,4)  = atan2(norm(cross(vector(i,:,2),vector(i,:,4))),dot(vector(i,:,2),vector(i,:,4)));
    theta(i,:,5)  = atan2(norm(cross(vector(i,:,2),vector(i,:,5))),dot(vector(i,:,2),vector(i,:,5)));
    theta(i,:,6)  = atan2(norm(cross(vector(i,:,4),vector(i,:,5))),dot(vector(i,:,4),vector(i,:,5)));
    theta(i,:,7)  = atan2(norm(cross(vector(i,:,1),vector(i,:,4))),dot(vector(i,:,1),vector(i,:,4)));
    theta(i,:,8)  = atan2(norm(cross(vector(i,:,4),vector(i,:,6))),dot(vector(i,:,4),vector(i,:,6)));
    theta(i,:,9)  = atan2(norm(cross(vector(i,:,1),vector(i,:,6))),dot(vector(i,:,1),vector(i,:,6)));
    theta(i,:,10) = atan2(norm(cross(vector(i,:,3),vector(i,:,6))),dot(vector(i,:,3),vector(i,:,6)));
    theta(i,:,11) = atan2(norm(cross(vector(i,:,3),vector(i,:,5))),dot(vector(i,:,3),vector(i,:,5)));
    theta(i,:,12) = atan2(norm(cross(vector(i,:,5),vector(i,:,6))),dot(vector(i,:,6),vector(i,:,6)));
end

% Identifiy the triangle for which the sum of the std of the related angles
% is minimal
for i = 1:size(theta,3)
    sdtheta(i)    = std(theta(:,:,i));
end
sdtriangle(1) = sum(sdtheta(1:3));
sdtriangle(2) = sum(sdtheta(4:6));
sdtriangle(3) = sum(sdtheta(7:9));
sdtriangle(4) = sum(sdtheta(10:12));
[~,ntriangle] = min(sdtriangle);

% Store the related marker and angle indices 
if ntriangle == 1
    nmarker = [1 2 3 4]; % Fourth marker listed only to recreate the last marker for inverse kinematics
    ntheta  = [1 2 3];
elseif ntriangle == 2
    nmarker = [1 3 4 2];
    ntheta  = [4 5 6];
elseif ntriangle == 3
    nmarker = [1 2 4 3];
    ntheta  = [7 8 9];
elseif ntriangle == 4
    nmarker = [2 3 4 1];
    ntheta  = [10 11 12];
end

% -------------------------------------------------------------------------
% STORE THE SOLID TRIANGLE (REFERENCE)
% -------------------------------------------------------------------------

% Compute the mean trajectory for each marker associated to the previously
% identified triangle
mtheta(1,1,1) = mean(theta(:,:,ntheta(1)),1);
mtheta(1,1,2) = mean(theta(:,:,ntheta(2)),1);
mtheta(1,1,3) = mean(theta(:,:,ntheta(3)),1);

% Keep only 75% of the coordinates for which the difference with the mean
% is minimal
nframe = 1:size(M1,1);
for i = 1:size(M1,1)
    test(i) = (mtheta(1,1,1)-theta(i,:,1))^2 + ...
              (mtheta(1,1,2)-theta(i,:,2))^2 + ...
              (mtheta(1,1,3)-theta(i,:,3))^2;
end
[~,ind] = max(test);
nframe = nframe(nframe~=ind);
clear ind;
% Compute the new mean with the "best" coordinates, for angles and markers
stheta(1,1,1)  = mean(theta(nframe,:,ntheta(1)),1);
stheta(1,1,2)  = mean(theta(nframe,:,ntheta(2)),1);
stheta(1,1,3)  = mean(theta(nframe,:,ntheta(3)),1);
smarker(1,:,1) = mean(marker(nframe,:,nmarker(1)),1);
smarker(1,:,2) = mean(marker(nframe,:,nmarker(2)),1);
smarker(1,:,3) = mean(marker(nframe,:,nmarker(3)),1);
smarker(1,:,4) = mean(marker(nframe,:,nmarker(4)),1); % Only for inverse kinematics

% Recreate the fourth marker (needed for inverse kinematics)

% Sort angles from the least deformed (lower std) to the most deformed
% (higher std) and sort markers accordingly
[stheta,ind]  = sort(stheta,3);
% nmarker       = [nmarker(ind(:,:,1)) nmarker(ind(:,:,2)) nmarker(ind(:,:,3))];
nmarker       = [nmarker(ind(:,:,1)) nmarker(ind(:,:,2)) nmarker(ind(:,:,3)) nmarker(4)];

% Define the vertices (ais) and normalised vertices (pis) of the solid 
% triangle (normalised = vertice at a unit distance from p1s)
a1s           = smarker(1,:,nmarker(1));
a2s           = smarker(1,:,nmarker(2));
a3s           = smarker(1,:,nmarker(3));
a4s           = smarker(1,:,nmarker(4)); % Only for inverse kinematics
p1s           = a1s;
p2s           = (a2s-a1s)/norm((a2s-a1s))+p1s;
p3s           = (a3s-a1s)/norm((a3s-a1s))+p1s;
clear ind;

% %% TEST
% figure;
% hold on;
% plot3(a1s(:,1),a1s(:,2),a1s(:,3),'Marker','x','Markersize',20,'Color','blue');
% plot3(a2s(:,1),a2s(:,2),a2s(:,3),'Marker','x','Markersize',20,'Color','blue');
% plot3(a3s(:,1),a3s(:,2),a3s(:,3),'Marker','x','Markersize',20,'Color','blue');
% plot3(p1s(:,1),p1s(:,2),p1s(:,3),'Marker','x','Markersize',20,'Color','red');
% plot3(p2s(:,1),p2s(:,2),p2s(:,3),'Marker','x','Markersize',20,'Color','red');
% plot3(p3s(:,1),p3s(:,2),p3s(:,3),'Marker','x','Markersize',20,'Color','red');
% line([a1s(:,1),a2s(:,1)],[a1s(:,2),a2s(:,2)],[a1s(:,3),a2s(:,3)],'Color','green');
% line([a1s(:,1),a3s(:,1)],[a1s(:,2),a3s(:,2)],[a1s(:,3),a3s(:,3)],'Color','green');
% %% TEST END

% -------------------------------------------------------------------------
% FIT THE SOLID TRIANGLE TO EACH MEASURED TRIANGLES
% -------------------------------------------------------------------------
for i = 1:size(M1,1)
%     disp(num2str(i*100/size(M1,1)));
    
    % Define the vertices (bi) and normalised vertices (qi) of the current
    % measured triangle
    b1(i,:)                    = marker(i,:,nmarker(1));
	b2(i,:)                    = marker(i,:,nmarker(2));
	b3(i,:)                    = marker(i,:,nmarker(3));
	b4(i,:)                    = marker(i,:,nmarker(4)); % Only for inverse kinematics
	q1(i,:)                    = b1(i,:);
	q2(i,:)                    = (b2(i,:)-b1(i,:))/norm((b2(i,:)-b1(i,:)))+q1(i,:);
	q3(i,:)                    = (b3(i,:)-b1(i,:))/norm((b3(i,:)-b1(i,:)))+q1(i,:);	
    
    % Compute the rigid transformation between the solid triangle and the
    % measured triangle by singular value decomposition algorithm
    % Söderkvist and Wedin (1993)
	[R(:,:,i),d(:,:,i),rms(i)] = soder([p1s;p2s;p3s],[q1(i,:);q2(i,:);q3(i,:)]);
    
    % Compute the solidified vertices (bis)
	b1s(i,:)                   = (R(:,:,i)*a1s'+d(:,:,i))';
    b2s(i,:)                   = (R(:,:,i)*a2s'+d(:,:,i))';
    b3s(i,:)                   = (R(:,:,i)*a3s'+d(:,:,i))';
    b4s(i,:)                   = (R(:,:,i)*a4s'+d(:,:,i))'; % Only for inverse kinematics
    
    % Update the marker coordinates
	eval(['M',num2str(nmarker(1)),'(i,:)=b1s(i,:);']);
	eval(['M',num2str(nmarker(2)),'(i,:)=b2s(i,:);']);
	eval(['M',num2str(nmarker(3)),'(i,:)=b3s(i,:);']);
	eval(['M',num2str(nmarker(4)),'(i,:)=b4s(i,:);']); % Only for inverse kinematics
    
    % Remove the cluster marker not used in the solid triangle
%     if isempty(find(nmarker==1))
%         M1 = NaN(size(M1));
%     elseif isempty(find(nmarker==2))
%         M2 = NaN(size(M2));
%     elseif isempty(find(nmarker==3))
%         M3 = NaN(size(M3));
%     elseif isempty(find(nmarker==4))
%         M4 = NaN(size(M4));
%     end
	
% 	%% TEST
% 	figure;
% 	hold on;
% 	plot3(b1(i,1),b1(i,2),b1(i,3),'Marker','x','Markersize',20,'Color','blue');
% 	plot3(b2(i,1),b2(i,2),b2(i,3),'Marker','x','Markersize',20,'Color','blue');
% 	plot3(b3(i,1),b3(i,2),b3(i,3),'Marker','x','Markersize',20,'Color','blue');
% 	plot3(q1(i,1),q1(i,2),q1(i,3),'Marker','x','Markersize',20,'Color','red');
% 	plot3(q2(i,1),q2(i,2),q2(i,3),'Marker','x','Markersize',20,'Color','red');
% 	plot3(q3(i,1),q3(i,2),q3(i,3),'Marker','x','Markersize',20,'Color','red');
% 	plot3(b1s(i,1),b1s(i,2),b1s(i,3),'Marker','x','Markersize',20,'Color','green');
% 	plot3(b2s(i,1),b2s(i,2),b2s(i,3),'Marker','x','Markersize',20,'Color','green');
% 	plot3(b3s(i,1),b3s(i,2),b3s(i,3),'Marker','x','Markersize',20,'Color','green');
% 	%% TEST END
end