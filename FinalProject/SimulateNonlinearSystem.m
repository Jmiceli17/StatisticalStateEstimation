% SimulatNonlinearSystem.m
%
% Description:
%    Simulate the nonlinear system and return time, states, and
%    measurements without using noise
%
% Inputs:
%   tspan - time span to conduct the simulation for, specified as
%           [0:step_size:final_time
%   x0 - initial condition
%   nGrndStations - number of ground stations to utilize
% Outputs:
%   t_vals - array of time values
%   x_vals - array of states obtained by solving the nonlinear differential
%            equations
%   y_vals - array of nonlinear measurements obtained using the state, note
%            that every time step has a corresponding measurement for each ground
%            station (some ground stations will produce NaN at a given time step)


function [t_vals, x_vals, y_vals] = SimulateNonlinearSystem(tspan, x_0, nGrndStations)


% Numerically integrate non-linear equations of motion over the provided
% time interval using the provided initial condition
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
u = [];
w = [];
[t_vals, x_vals] = ode45(@(t,y) two_body_diff_eq(t, y, u, w), tspan, x_0, opts);
% [t_vals, x_vals] = ode45(@(t,y) nonLinearOde(t, y, u, w), tspan, x_0, opts);


% Generate measurements for each ground station
y_vals = zeros(size(t_vals,1),3,nGrndStations);
visible_gs_array = zeros(size(t_vals,1),nGrndStations);
gsState_vals = zeros(size(t_vals,1), 4, nGrndStations);

% Loop over each ground station
for gsIdx = 1:nGrndStations
    
    for t = 1:size(t_vals,1)-1
        % Get the time and the state
        time = t_vals(t);
        X = x_vals(t,:)';
        
        % Get Ground station state
        gsState = GetGroundStationState(time, gsIdx);
        
        % Compute the measurement using full non-linear equations
        evaluate_visibility = true;
        Y = nonLinearMeasurementOde(time, X, gsState, evaluate_visibility);
      
        % Log the ground station state
        gsState_vals(t, :, gsIdx) = gsState;
        
        % Log the measurement
        y_vals(t+1, :, gsIdx) = Y';
        
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
    


% ---- Actual State Plot ----
figure()
suptitle('States vs Time, Nonlinear Dynamics Simulation');
set(findall(gcf,'type','text'),'FontSize',16)


subplot(4,1,1);
p1 = plot(t_vals, x_vals(:,1), 'b-', 'linewidth',  2);
xlabel('Time (s)','FontSize', 12)
ylabel('x_{1}(t) [km]','FontSize', 12)
grid on
subplot(4,1,2)
p2 = plot(t_vals, x_vals(:,2), 'r-', 'linewidth', 2);
xlabel('Time (s)','FontSize', 12)
ylabel('x_{2}(t) [km/s]','FontSize', 12)
grid on
subplot(4,1,3)
p3 = plot(t_vals, x_vals(:,3), 'g-', 'linewidth', 2);
xlabel('Time (s)','FontSize', 12)
ylabel('x_{3}(t) [km]','FontSize', 12)
grid on
subplot(4,1,4)
p4 = plot(t_vals, x_vals(:,4), 'k-', 'linewidth', 2);
xlabel('Time (s)','FontSize', 12)
ylabel('x_{4}(t) [km/s]','FontSize', 12)
grid on
set(findall(gcf,'type','line'),'linewidth',2)

% % ---- Plot the ground station states ----
% figure()
% suptitle('Nonlinear Simulation - States of Ground Stations in XY Reference Frame');
% set(findall(gcf,'type','text'),'FontSize',16)
% for gsIdx = 1:nGrndStations
%     
%     % DEBUGGING (only plot first gs)
%     if gsIdx == 1
%         
%         % Generate random color scheme to use for this ground station's data
%         rand_color = [rand, rand, rand];
% 
%         subplot(4,1,1);
%         hold on;
%         plot(t_vals, gsState_vals(:,1,gsIdx), 'x', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('X_{s} [km]','FontSize', 12)
%         grid on
% 
%         subplot(4,1,2)
%         hold on;
%         plot(t_vals, gsState_vals(:,2,gsIdx), 'x', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('$\dot{X_{s}} [km/s]$','FontSize', 12, 'Interpreter', 'latex')
%         grid on
% 
%         subplot(4,1,3)
%         hold on;
%         plot(t_vals, gsState_vals(:,3,gsIdx), 'x', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('Y_{s} [km]','FontSize', 12)
%         grid on
% 
%         subplot(4,1,4)
%         hold on;
%         plot(t_vals, gsState_vals(:,4,gsIdx), 'X', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('$\dot{Y_{s}} [km/s]$','FontSize', 12, 'Interpreter', 'latex')
%         grid on
%     end
% end


% ---- Measurement State Plot ----
figure()
suptitle('Full Nonlinear Model Data Simulation');
set(findall(gcf,'type','text'),'FontSize',18)


for gsIdx = 1:nGrndStations
    
    % Generate random color scheme to use for this ground station's data
    rand_color = [rand, rand, rand];
    
    subplot(4,1,1);
    hold on;
    plot(t_vals, y_vals(:,1,gsIdx), 'x', 'Color', rand_color);
    hold off;
    xlabel('Time (s)','FontSize', 12)
    ylabel('$\rho(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
    grid on
    
    subplot(4,1,2)
    hold on;
    plot(t_vals, y_vals(:,2,gsIdx), 'x', 'Color', rand_color);
    hold off;
    xlabel('Time (s)','FontSize', 12)
    ylabel('$\dot{\rho}(t) [km/s]$','FontSize', 14, 'Interpreter', 'latex')
    grid on
    
    subplot(4,1,3)
    hold on;
    plot(t_vals, y_vals(:,3,gsIdx), 'x', 'Color', rand_color);
    hold off;
    xlabel('Time (s)','FontSize', 12)
    ylabel('$\phi(t) [rad]$','FontSize', 14, 'Interpreter', 'latex')
    grid on
    
    subplot(4,1,4)
    hold on;
    plot(t_vals, visible_gs_array(:,gsIdx), '*', 'Color', rand_color);
    hold off;
    xlabel('Time (s)','FontSize', 12)
    ylabel('Visible Ground Station ID','FontSize', 12)
    grid on
    
end

end