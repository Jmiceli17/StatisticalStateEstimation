% SimulateLinearizedSystem.m
% 
% Description:
%   Simulates the linearized system of equations for the the 2D satellite
%   orbit determination problem and plots the results.
% 
% Inputs:
%   A_nom     - Jacobian of nonlinear system wrt state (dF/dx) evaluated at
%               nominal state (in this case, nominal state is function of time)
%   B_nom     - Jacobian of nonlinear system wrt controls (dF/dU),
%               time-invariant in this case
%   x_pert_0  - initial perturbed state
%   u_pert_0  - initial perturbed inputs
%   num_steps - number of time steps to simulate
% 
% Outputs: 
%   time_vals - double array of time values 
%   x         - array of estimated full state
%   x_nom_vals - array of nominal state
%   y_nom_vals - array of linearized measurements of nominal state


function [time_vals, x, xnom_evals, y_nom_vals] = SimulateLinearizedSystem(x_pert_0, u_pert_0, num_steps, num_gs)

global xinit num_states num_inputs delta_t  

% Initialize simulation time
time_vals = zeros(num_steps+1,1);

% Initialize perturbation state, controls
x_pert = zeros(num_steps+1, num_states);
x_pert(1,:) = x_pert_0';

u_pert = zeros(num_steps+1, num_inputs);
u_pert(1,:) = u_pert_0';

% Initialize Nominal State
xnom_evals = zeros(num_steps+1, num_states);
% xnom_evals(1,:) = xinit;
xnom_evals(1,:) = GetNominalState(0.0)';

% Initialize actual state
x = zeros(num_steps+1, num_states);
x(1,:) = xinit + x_pert_0;

for k = 1:num_steps
    % Current simulation time
    sim_time = double(k*delta_t);
    
    % Get the curernt nominal state 
    xnom_k = xnom_evals(k,:)';
    
    % Evaluate symbolic jacobian matrices
    unom = zeros(2,1);      % zeros for now
    wprocess = zeros(4,1);  % zeros for now
    [A_nom_eval, B_nom_eval] = LinearizedDynamicsOde(xnom_k, unom, wprocess);
    
    % Compute the linearized DT system matrices
    F_tilde_k = eye(num_states, num_states) + delta_t * A_nom_eval;
    G_tilde_k = delta_t * B_nom_eval;
    
    % Previous perturbation state and inputs
    % Note k = 1 corresponds to the initial condtion
    x_pert_k = x_pert(k,:)';    
    u_pert_k = u_pert(k,:)';
    
    % Compute new perturbation state using linearized sys parameters
    x_pert_kplus1 = F_tilde_k*x_pert_k + G_tilde_k*u_pert_k;

    % Compute the new nominal state (circular orbit with constant ang vel)
    xnom_kplus1 = GetNominalState(sim_time);
    
    % Compute new total state
    x_kplus1 = x_pert_kplus1 + xnom_kplus1; 
    
    % Log new data
    time_vals(k+1)    = sim_time;
    xnom_evals(k+1,:) = xnom_kplus1';
    x_pert(k+1,:)     = x_pert_kplus1';
    x(k+1,:)          = x_kplus1';
end


% Arrays to hold perturbation, nominal, and total measurement vectors
y_pert_vals      = zeros(num_steps+1, 3, num_gs);
y_nom_vals       = zeros(num_steps+1, 3, num_gs);
y_vals           = zeros(num_steps+1, 3, num_gs);
visible_gs_array = zeros(num_steps+1, num_gs);
gsState_vals     = zeros(num_steps+1, 4, num_gs);




% Loop over each ground station to get the linearized measurements
for gsIdx = 1:num_gs
    
    % Initialize ground station state array for this ground station
    gsState_vals(1, :, gsIdx) = GetGroundStationState(0.0, gsIdx);

    
    for k = 1:num_steps
              
        sim_time      = (k+1)*delta_t;
        x_nom_kplus1  = xnom_evals(k+1,:)';  % nominal state
        x_pert_kplus1 = x_pert(k+1,:)';      % perturbed state
        u_unom_kplus1 = zeros(2,1);          % zeros for now
        v_kplus1      = zeros(3,1);          % zeros for now
        x_total       = x(k+1,:)';           % total state
        
        % Get the state of the ground station
        gsState = GetGroundStationState(sim_time, gsIdx);
        
        elev_angle = atan2( (x_total(3)-gsState(3)),(x_total(1) - gsState(1)) );
        theta = atan2(gsState(3),gsState(1)); % returns values [-pi,pi], atan returns [-pi/2,pi/2]
        err = min(2*pi - abs(elev_angle - theta), abs(elev_angle - theta) );
    
        if err <= pi/2
            % Compute the measurement using linearized sensing matrix (note,
            % the linearized CT sensing matrix is equivalent to the DT sensing
            % matrix, H_tilde = C_tilde)
            H_tilde_kplus1 = LinearizedMeasurementOde(x_nom_kplus1, u_unom_kplus1, v_kplus1, gsState);
            y_pert_kplus1  = H_tilde_kplus1 * x_pert_kplus1;    % Perturbed meas
            
            % Compute the nominal measurement
            evaluate_visibility = true;
            y_nom_kplus1   = nonLinearMeasurementOde(sim_time, x_nom_kplus1, gsState, evaluate_visibility);
            
            % Wrap elevation angle to -pi to pi
            y_pert_kplus1(3) = wrapToPi(y_pert_kplus1(3));
            y_nom_kplus1(3)  = wrapToPi(y_nom_kplus1(3));
            
            y_kplus1       = y_pert_kplus1 + y_nom_kplus1;
            
            % Wrap elevation angle to -pi to pi
            y_kplus1(3)  = wrapToPi(y_kplus1(3));
            
            
            % Log the measurement
            y_pert_vals(k+1, :, gsIdx) = y_pert_kplus1';
            y_nom_vals(k+1, :, gsIdx)  = y_nom_kplus1';
            y_vals(k+1, :, gsIdx)      = y_kplus1';

        else
            y_pert_kplus1  = [NaN, NaN, NaN];   
            y_nom_kplus1   = [NaN, NaN, NaN];
            y_kplus1 = [NaN, NaN, NaN];
            
            y_pert_vals(k+1, :, gsIdx) = y_pert_kplus1';
            y_nom_vals(k+1, :, gsIdx)  = y_nom_kplus1';
            y_vals(k+1, :, gsIdx)      = y_kplus1';
        end
        
        % Track when the GS is visible for plotting
        if isnan(y_pert_kplus1)
            visible_gs_array(k+1,gsIdx) = NaN;
        else
            visible_gs_array(k+1,gsIdx) = gsIdx;
        end
        
        
        % Log the state of the ground station
        gsState_vals(k+1, :, gsIdx) = gsState';
        
    end
end



% % % % ---- Actual and Nominal State Plot ----
% % % figure(1)
% % % suptitle('States vs Time, Linearized Approximate DT Dynamics Simulation');
% % % set(findall(gcf,'type','text'),'FontSize',18)
% % % 
% % % 
% % % subplot(4,1,1);
% % % p1 = plot(time_vals, x(:,1), 'b-', 'linewidth',  2);
% % % hold on;
% % % plot(time_vals, xnom_evals(:,1), 'b--', 'linewidth',  2);
% % % legend('Actual', 'Nominal')
% % % xlabel('Time (s)','FontSize', 14)
% % % ylabel('$x(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% % % grid on
% % % subplot(4,1,2)
% % % p2 = plot(time_vals, x(:,2), 'r-', 'linewidth', 2);
% % % hold on;
% % % plot(time_vals, xnom_evals(:,2), 'r--', 'linewidth',  2);
% % % legend('Actual', 'Nominal')
% % % xlabel('Time (s)','FontSize', 14)
% % % ylabel('$\dot{x}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % % grid on
% % % 
% % % subplot(4,1,3)
% % % p3 = plot(time_vals, x(:,3), 'g-', 'linewidth', 2);
% % % hold on;
% % % plot(time_vals, xnom_evals(:,3), 'g--', 'linewidth',  2);
% % % legend('Actual', 'Nominal')
% % % xlabel('Time (s)','FontSize', 14)
% % % ylabel('$y(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% % % grid on
% % % 
% % % subplot(4,1,4)
% % % p4 = plot(time_vals, x(:,4), 'k-', 'linewidth', 2);
% % % hold on;
% % % plot(time_vals, xnom_evals(:,4), 'k--', 'linewidth',  2);
% % % legend('Actual', 'Nominal')
% % % xlabel('Time (s)','FontSize', 14)
% % % ylabel('$\dot{y}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % % grid on
% % % set(findall(gcf,'type','line'),'linewidth',2)
% % % 
% % % % ---- Perturbation State Plot ----
% % % figure(2)
% % % suptitle('Linearized Approximation Perturbation State vs Time');
% % % set(findall(gcf,'type','text'),'FontSize',18)
% % % 
% % % 
% % % subplot(4,1,1);
% % % p1 = plot(time_vals, x_pert(:,1), 'b-', 'linewidth',  2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('$\delta x(t) [km]$','FontSize', 14,'Interpreter', 'latex')
% % % grid on
% % % subplot(4,1,2)
% % % p2 = plot(time_vals, x_pert(:,2), 'r-', 'linewidth', 2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('$\delta\dot{x}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % % grid on
% % % subplot(4,1,3)
% % % p3 = plot(time_vals, x_pert(:,3), 'g-', 'linewidth', 2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('$\delta y(t) [km]$','FontSize', 14,'Interpreter', 'latex')
% % % grid on
% % % subplot(4,1,4)
% % % p4 = plot(time_vals, x_pert(:,4), 'k-', 'linewidth', 2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('$\delta\dot{y}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % % grid on
% % % set(findall(gcf,'type','line'),'linewidth',2)
% % % 
% % % 
% % % % ---- Nominal State Plot ----
% % % figure(3)
% % % suptitle('Nominal State vs Time');
% % % set(findall(gcf,'type','text'),'FontSize',18)
% % % 
% % % 
% % % subplot(4,1,1);
% % % p1 = plot(time_vals, xnom_evals(:,1), 'b-', 'linewidth',  2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('x_{1}(t) [km]','FontSize', 12)
% % % grid on
% % % subplot(4,1,2)
% % % p2 = plot(time_vals, xnom_evals(:,2), 'r-', 'linewidth', 2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('x_{2}(t) [km/s]','FontSize', 12)
% % % grid on
% % % subplot(4,1,3)
% % % p3 = plot(time_vals, xnom_evals(:,3), 'g-', 'linewidth', 2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('x_{3}(t) [km]','FontSize', 12)
% % % grid on
% % % subplot(4,1,4)
% % % p4 = plot(time_vals, xnom_evals(:,4), 'k-', 'linewidth', 2);
% % % xlabel('Time (s)','FontSize', 12)
% % % ylabel('x_{4}(t) [km/s]','FontSize', 12)
% % % grid on
% % % set(findall(gcf,'type','line'),'linewidth',2)

% % ---- Plot the ground station states ----
% figure()
% suptitle('Linearized Simulation - States of Ground Stations in XY Reference Frame');
% set(findall(gcf,'type','text'),'FontSize',16)
% for gsIdx = 1:num_gs
%     
%     % DEBUGGING (only plot first gs)
%     if gsIdx == 1
%         
%         % Generate random color scheme to use for this ground station's data
%         rand_color = [rand, rand, rand];
% 
%         subplot(4,1,1);
%         hold on;
%         plot(time_vals, gsState_vals(:,1,gsIdx), 'x', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('X_{s} [km]','FontSize', 12)
%         grid on
% 
%         subplot(4,1,2)
%         hold on;
%         plot(time_vals, gsState_vals(:,2,gsIdx), 'x', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('$\dot{X_{s}} [km/s]$','FontSize', 12, 'Interpreter', 'latex')
%         grid on
% 
%         subplot(4,1,3)
%         hold on;
%         plot(time_vals, gsState_vals(:,3,gsIdx), 'x', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('Y_{s} [km]','FontSize', 12)
%         grid on
% 
%         subplot(4,1,4)
%         hold on;
%         plot(time_vals, gsState_vals(:,4,gsIdx), 'X', 'Color', rand_color);
%         hold off;
%         xlabel('Time (s)','FontSize', 12)
%         ylabel('$\dot{Y_{s}} [km/s]$','FontSize', 12, 'Interpreter', 'latex')
%         grid on
%     end
% end


% % % % ---- Measurement State Plot ----
% % % figure()
% % % suptitle('Linearized Model Data Simulation');
% % % set(findall(gcf,'type','text'),'FontSize',18)
% % % 
% % % 
% % % for gsIdx = 1:num_gs
% % %     
% % %     % Generate random color scheme to use for this ground station's data
% % %     rand_color = [rand, rand, rand];
% % %     
% % %     subplot(4,1,1);
% % %     hold on;
% % %     plot(time_vals, y_vals(:,1,gsIdx), 'x', 'Color', rand_color);
% % %     hold off;
% % %     xlabel('Time (s)','FontSize', 12)
% % %     ylabel('$\rho(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% % %     grid on
% % %     
% % %     subplot(4,1,2)
% % %     hold on;
% % %     plot(time_vals, y_vals(:,2,gsIdx), 'x', 'Color', rand_color);
% % %     hold off;
% % %     xlabel('Time (s)','FontSize', 12)
% % %     ylabel('$\dot{\rho}(t) [km/s]$','FontSize', 14, 'Interpreter', 'latex')
% % %     grid on
% % %     
% % %     subplot(4,1,3)
% % %     hold on;
% % %     plot(time_vals, y_vals(:,3,gsIdx), 'x', 'Color', rand_color);
% % %     hold off;
% % %     xlabel('Time (s)','FontSize', 12)
% % %     ylabel('$\phi(t) [rad]$','FontSize', 14, 'Interpreter', 'latex')
% % %     grid on
% % %     
% % %     subplot(4,1,4)
% % %     hold on;
% % %     plot(time_vals, visible_gs_array(:,gsIdx), '*', 'Color', rand_color);
% % %     hold off;
% % %     xlabel('Time (s)','FontSize', 12)
% % %     ylabel('Visible Ground Station ID','FontSize', 14)
% % %     grid on
% % %     
% % % end

end