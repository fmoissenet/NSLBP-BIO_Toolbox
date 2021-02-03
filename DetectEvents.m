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

function Trial = DetectEvents(Trial,vec1,vec2,type,threshold)

% Type 1: TO BE DEFINED
if type == 1
    for time = 1:size(vec1,1)
        temp = atan2(norm(cross(vec1(time,:),vec2(time,:))),dot(vec1(time,:),vec2(time,:))); % Angle in radians
        if rad2deg(temp) < threshold
            a(time) = temp;
        else
            a(time) = NaN;
        end
        clear temp;
    end
    temp1 = [];
    temp2 = [];
    for time = 1:size(a,2)-1
        if ~isnan(a(time)) && isnan(a(time+1))
            temp1 = [temp1 time];
        elseif isnan(a(time)) && ~isnan(a(time+1))
            temp2 = [temp2 time];
        end
    end
    ind1 = [];
    ind2 = [];
    ind1 = [ind1 10];
    for j = 1:size(temp1,2)
        ind2 = [ind2 temp1(j)];
        if j < size(temp1,2)
            ind1 = [ind1 fix((temp2(j)+temp2(j+1))/2)];              
        end
    end
    ind1 = [ind1 Trial.n1-10];
    Trial.Event(1).label = 'start';
    Trial.Event(1).value = ind1(1:size(ind1,2));
    Trial.Event(2).label = 'back';
    Trial.Event(2).value = ind2(1:size(ind2,2));
    clear vec1 vec2 time a ind1 ind2 temp1 temp2;

% Type 2: The back event appears when the angle used for detection is close
% to its maximum value, and the start value to its minimum value   
elseif type == 2
    for time = 1:size(vec1,1)
        temp = atan2(norm(cross(vec1(time,:),vec2(time,:))),dot(vec1(time,:),vec2(time,:))); % Angle in radians
        if rad2deg(temp) < threshold
            a(time) = temp;
        else
            a(time) = NaN;
        end
        clear temp;
    end
    temp1 = [];
    temp2 = [];
    for time = 1:size(a,2)-1
        if ~isnan(a(time)) && isnan(a(time+1))
            temp1 = [temp1 time];
        elseif isnan(a(time)) && ~isnan(a(time+1))
            temp2 = [temp2 time];
        end
    end
    ind1 = [];
    ind2 = [];
    ind1 = [ind1 10];
    for j = 1:size(temp1,2)
        if j < size(temp1,2)
            ind1 = [ind1 fix((temp2(j)+temp1(j+1))/2)];
        end
        ind2 = [ind2 fix((temp1(j)+temp2(j))/2)];
    end
    ind1 = [ind1 Trial.n1-10];
    Trial.Event(1).label = 'start';
    Trial.Event(1).value = ind1(1:size(ind1,2));
    Trial.Event(2).label = 'back';
    Trial.Event(2).value = ind2(1:size(ind2,2));
    clear vec1 vec2 time a ind1 ind2 temp1 temp2;

% Type 3: The angle used for detection reach its maximum value between
% start and back
elseif type == 3
    for time = 1:size(vec1,1)
        temp = atan2(norm(cross(vec1(time,:),vec2(time,:))),dot(vec1(time,:),vec2(time,:))); % Angle in radians
        if rad2deg(temp) < threshold
            a(time) = temp;
        else
            a(time) = NaN;
        end
        clear temp;
    end
    temp1 = [];
    temp2 = [];
    for time = 1:size(a,2)-1
        if ~isnan(a(time)) && isnan(a(time+1))
            temp1 = [temp1 time];
        elseif isnan(a(time)) && ~isnan(a(time+1))
            temp2 = [temp2 time];
        end
    end
    ind1 = [];
    ind2 = [];
    ind1 = [ind1 10];
    for j = 1:size(temp1,2)
        if j < size(temp1,2)
            ind1 = [ind1 fix((temp2(j)+temp1(j+1))/2)];
        end
        ind2 = [ind2 fix((temp1(j)+temp2(j))/2)];
    end
    ind1 = [ind1 Trial.n1-10];
    Trial.Event(1).label = 'start';
    Trial.Event(1).value = ind1(1:2:size(ind1,2));
    Trial.Event(2).label = 'back';
    Trial.Event(2).value = ind1(2:2:size(ind1,2));
    clear vec1 vec2 time a ind1 ind2 temp1 temp2;
end