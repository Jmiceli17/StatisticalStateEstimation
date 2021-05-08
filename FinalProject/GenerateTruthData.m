% GenerateTruthData.m
%
% Description:
%   Generates noisy measurements, states, and nominal measurements 
% Inputs:
%   useNoise - boolean indicating if noise should be used or not
%   num_steps - int numbe of time steps to generate truth data for
%   nGrndStations - number of ground stations to evaluate
%   x0 - true initial state of the s/c (x,xdot, y, ydot)
%   Q_truth - process noise intensity covariance
%   R_truth - measurement noise intensity covariance
%
% Outputs:
%   t_vals - array of time values [num_steps + 1, 1]
%   x_truth_vals - Array of true noisy state of the spacecraft for 0 to 
%                  num_steps time steps
%   y_truth_vals - Array of measurements of the true state of the s/c
%                  produced by each ground station (i.e. at every time 
%                  step, there is a measurement from all ground stations, 
%                  could be a 3x1 array or 3x1 NaN
%   array)
%   y_nom_truth_vals - Array of meaasurements of the NOMINAL state of the
%                      s/c
%   visible_gs_array - Array indicating which ground stations are visible
%                      at this time step


function [t_vals, x_truth_vals, y_truth_vals, y_nom_truth_vals, visible_gs_array] = GenerateTruthData(useNoise, num_steps, nGrndStations, x0, Q_truth, R_truth)

    global delta_t num_states

    % ODE45 options
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % No inputs
    u = [];

    t_vals            = zeros(num_steps + 1,1);
    x_truth_vals      = zeros(num_steps + 1, num_states);
    x_truth_vals(1,:) = x0;
    
    % Generate noise values
    Sw = chol(Q_truth, 'lower');
    wtilde = useNoise * (Sw*(randn(length(t_vals),2))')';
    
    % Break up the integration of the nonlinear dynamics so different noise
    % can be applied at each time step
    for k = 1:num_steps
        
        % Sample the noise vector
        wk = wtilde(k,:)';

        % Integrate the spacecraft dynamics with noise added as ZOH input
%         [T_kplus1, x_kplus1] = ode45(@(t,y) two_body_diff_eq(t, y, u, wk), [0 delta_t], x0, opts);
        [~, x_kplus1] = ode45(@(t,y) nonLinearOde(t, y, u, wk), [0:1:delta_t], x0, opts);
        
        % Update inititial state to the final state of the solution
        x0 = x_kplus1(end,:)'; 
        
        % Store the results
        x_truth_vals(k+1,:) =  x_kplus1(end,:)';
        t_vals(k+1) = k*delta_t;
    end
   
    
    % Generate measurements for each ground station
    y_truth_vals     = zeros(size(t_vals,1),3,nGrndStations);
    y_nom_truth_vals = zeros(size(t_vals,1),3,nGrndStations);
    visible_gs_array = zeros(size(t_vals,1),nGrndStations);
    gsState_vals     = zeros(size(t_vals,1), 4, nGrndStations);

    % Loop over each ground station
    for gsIdx = 1:nGrndStations

        for t = 1:size(t_vals,1)-1
            % Get the time and the state
            time = t_vals(t);
            X = x_truth_vals(t,:)';
            xnom = GetNominalState(time);

            % Get Ground station state
            gsState = GetGroundStationState(time, gsIdx);

            % Compute noisey measurements using full non-linear equations
            evaluate_visibility = true;
            Y = nonLinearMeasurementOde(time, X, gsState, evaluate_visibility);
            ynom = nonLinearMeasurementOde(time, xnom, gsState, evaluate_visibility);

            % Generate noise value to apply to measurement
            Sv = chol(R_truth, 'lower');
            v = useNoise* (Sv*randn(3,1));

            % Log the ground station state
            gsState_vals(t, :, gsIdx) = gsState;

            % Log the measurement
            y_truth_vals(t+1, :, gsIdx) = (Y + v)';
            y_nom_truth_vals(t+1, :, gsIdx) = (ynom + v)';

            % Track when the GS is visible for plotting
            if isnan(Y)
                visible_gs_array(t,gsIdx) = NaN;
            else
                visible_gs_array(t,gsIdx) = gsIdx;
            end

        end

        % DEBUGGING
        % Log the last ground station state
        gsState = GetGroundStationState(t_vals(end), gsIdx);

        gsState_vals(end, :, gsIdx) = gsState;

    end
    
    
end

