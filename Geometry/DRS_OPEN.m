function [x, y] = DRS_OPEN(x,y,DRS_angle)

% DRS_OPEN Rotates the flap profile for the DRS ON configuration.
%
% Input:
%
%   x --> x-coordinates of the 2 profiles (by column)
%   y --> y-coordinates of the 2 profiles (by row)
%
% Example:
%
%   [x,y] = DRS_OPEN(x,y,alpha(2) - alpha(1))
%

    %% Unpack input vectors
    if ~iscell(x)
        
        %Main Profile Vectors
        x1 = x(:,1);
        y1 = y(:,1);
    
        %Flap Profile Vectors
        x2 = x(:,2);
        y2 = y(:,2);
    
    
    
        flap_DRS_OFF = [x2'; y2'];
    
        %% Rotation Matrix
        R = [cosd(DRS_angle) -sind(DRS_angle);...
             sind(DRS_angle)  cosd(DRS_angle)];
    
        %% Shift the center of roation to the hinge point (L.E.) 
    
        %New center definition
        x_center = x2(end);
        y_center = y2(end);
        center = repmat([x_center; y_center],1,length(x));
    
        %Shift reference to the hinge point
        flap_DRS_OFF = flap_DRS_OFF - center;
    
        %% Rotation along the hinge
        flap_DRS_ON = R * flap_DRS_OFF;
    
        flap_DRS_ON = flap_DRS_ON + center;
    
    
        %%Pack output vectors
        x2 = flap_DRS_ON(1,:)';
        y2 = flap_DRS_ON(2,:)';
    
        x = [x1, x2];
        y = [y1, y2];
    
    else
    %% Unpack Input Arrays
    
        x1 = x{1};
        y1 = y{1};
        
        x2 = x{2};
        y2 = y{2};
        
        flap_DRS_OFF = [x2'; y2'];
    
        %% Rotation Matrix
        R = [cosd(DRS_angle) -sind(DRS_angle);...
             sind(DRS_angle)  cosd(DRS_angle)];
    
        %% Shift the center of roation to the hinge point (L.E.) 
    
        %New center definition
        x_center = x2(end);
        y_center = y2(end);
        center = repmat([x_center; y_center],1,length(x2));
    
        %Shift reference to the hinge point
        flap_DRS_OFF = flap_DRS_OFF - center;
    
        %% Rotation along the hinge
        flap_DRS_ON = R * flap_DRS_OFF;
    
        flap_DRS_ON = flap_DRS_ON + center;
    
    
        %%Pack output vectors
        x2 = flap_DRS_ON(1,:)';
        y2 = flap_DRS_ON(2,:)';
    
        x = {x1, x2};
        y = {y1, y2};
    end

end