% Get the state of this ground station at this time

function gsState = GetGroundStationState(time, groundStationIdx)

    global re we


    % Position and velocity of ground station
    theta_0 = (groundStationIdx - 1)* pi/6;
    xs = re*cos(we*time + theta_0);
    xs_dot = -re*we*sin(we*time + theta_0);
    ys = re*sin(we*time + theta_0);
    ys_dot = re*we*cos(we*time + theta_0);
    gsState = [xs, xs_dot, ys, ys_dot]';    % State vector of the ground station 
    
    
    
    
end
