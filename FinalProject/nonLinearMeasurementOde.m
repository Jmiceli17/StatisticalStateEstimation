% nonLinearMeasurementOde.m 
%
% Description:
%   Nonlinear system of measurement equations assuming no sensing noise
%   this function is equivalent to the "H" function in standard nonlinear 
%   SS form

function Y = nonLinearMeasurementOde(~, X, gsState, evaluate_visibility )

    % Position and velocity of ground station    
    xs = gsState(1);
    xs_dot = gsState(2);
    ys = gsState(3);
    ys_dot = gsState(4);
    theta = atan2(ys,xs); % returns values [-pi,pi], atan returns [-pi/2,pi/2]

    
    
    % Build the outputs
    range = sqrt((X(1) - xs)^2 + (X(3) - ys)^2);
    range_rate = ((X(1) - xs)*(X(2) - xs_dot)+(X(3) - ys)*(X(4) - ys_dot))/range;
    elev_angle = atan2((X(3)-ys),(X(1) - xs));
    
    % If we're evaluating visibility of the gs
    if evaluate_visibility            
        % Check if GS is visible, if err is the absolute value of the angle
        % measuring from directly overhead to the satellite
        err = min(2*pi - abs(elev_angle - theta), abs(elev_angle - theta) );
    
        if err <= pi/2
            Y = [range;
            range_rate;
            elev_angle];
        else 
            Y = [NaN;
                NaN;
                NaN];
        end
        
    else
        Y = [range;
            range_rate;
            elev_angle];
        
    end

    

end