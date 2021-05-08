% RunExtendedKalmanFilter.m
% 
% Description:
%   Original EKF, designed to ingest "truth" data and report estimates, note
%   this filter is not set up to utilize the provided measurement data.

function [x_estim_vals, P_vals, NEES_vals, NIS_vals, estim_error_vals] = ...
    RunExtendedKalmanFilter(y_Truth, Q_true, R_true, x_Truth, xinit, P_0, num_steps)

    global num_states delta_t

    % ODE45 options and initialization
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    x_0_ode45 = xinit;  

    
    % Initialize arrays to store estimates and est err covariances
    x_estim_vals = zeros(num_steps+1, num_states);
    x_estim_vals(1, :) = xinit';
    
    t_vals = zeros(num_steps+1, 1);

    P_vals        = zeros(num_states, num_states, num_steps+1);
    P_vals(:,:,1) = P_0;
    
    estim_error_vals = zeros(num_steps+1, 4);
    
    NEES_vals    = zeros(num_steps+1,1);
    est_err_0    = (x_Truth(1, :)' - xinit);        % column vector
    NEES_vals(1) = est_err_0'*inv(P_0)*est_err_0;
    
    NIS_vals    = zeros(num_steps+1, 1);
    NIS_vals(1) = NaN;     % TODO: should this be NaN?  
    
    % Tuning parameters
    % HIGH Q -> low confidence in dynamics (trust the measurements)
    % HIGH R -> low confidence in sensors (trust the prediction step)
    %%Case 1: KF uses the correct/true Q and R values for DT system:
    Q_k = Q_true;
    R_kplus1 = R_true;
    
    for k =1:num_steps

        
        % Get values from t = tk
        x_k          = x_estim_vals(k,:)';
        P_k          = P_vals(:,:,k);
        u = zeros(2,1);
        w = zeros(2,1);
        v = zeros(3,1);
        
        % Linearized system parameters (linearized about the current
        % estimate)
        [A_tilde_k, ~] = LinearizedDynamicsOde(x_k, u, w);
        Omga_tilde_k = delta_t * [0 0; 1 0; 0 0; 0 1];  % Time invariant, Q = dT * Gamma(t)

        F_tilde_k = eye(num_states, num_states) + delta_t * A_tilde_k; 
        
        
        % --- Time update prediction step ---
        % Integrate nonlinea dynamics from tk to tk+1
       [~, x_kplus1_min] = ode45(@(t,y) nonLinearOde(t, y, u, w), [0:1:delta_t], x_0_ode45, opts);  
        x_0_ode45 = x_kplus1_min(end,:)';
        
        % Updated perturbation estimate and est error covariance
        x_kplus1_min = x_kplus1_min(end,:)'; % last element from ode45 integration solution
        
        P_kplus1_min  = F_tilde_k*P_k*F_tilde_k' + Omga_tilde_k*Q_k*Omga_tilde_k';

        % Get values from t = tk+1
        t_kplus1 = delta_t*(k+1);
        
        % Get the first available measurement at t=tk+1
        for gsIdx = 1:12
            % Get Ground station state
            gsState = GetGroundStationState(t_kplus1, gsIdx);
            % Get nonlinear measurement estimate
            evaluate_visibility = false;
            y_kplus1_min = nonLinearMeasurementOde(t_kplus1, x_kplus1_min, gsState, evaluate_visibility );
            y_tru_kplus1 = y_Truth(k+1, :, gsIdx)';
            
            % If measurement is not valid, skip
            if any(isnan(y_tru_kplus1))
            % skip
                x_kplus1  = x_kplus1_min;
                P_kplus1  = P_kplus1_min;


                NEES_vals(k+1) = NaN;
                NIS_vals(k+1) = NaN;
                
             
            else
                % Compute gain
                H_tilde_kplus1 = LinearizedMeasurementOde(x_kplus1_min, u, v, gsState);
                K_kplus1       = P_kplus1_min*H_tilde_kplus1'*inv(H_tilde_kplus1*P_kplus1_min*H_tilde_kplus1' + R_kplus1);
                innov_kplus1   = y_tru_kplus1 - y_kplus1_min;
                % Wrap to pi
                innov_kplus1(end) = atan2(sin(innov_kplus1(end)), cos(innov_kplus1(end)));
                
                % Updated total state estimate and estimate error covariance
                x_kplus1 = x_kplus1_min + K_kplus1*innov_kplus1;
                P_kplus1 = (eye(num_states, num_states) - K_kplus1*H_tilde_kplus1)*P_kplus1_min;

                % Compute and Store NEES and NIS Statistics
                est_err_kplus1 = (x_Truth(k+1,:)' - x_kplus1);  % column vector
                NEES_vals(k+1) = est_err_kplus1'*inv(P_kplus1)*est_err_kplus1;

                S_kplus1 = H_tilde_kplus1*P_kplus1_min*H_tilde_kplus1' + R_kplus1;

                % This guarantees that S is positive semi-definite (which it should be because it's a covariance matrix)
                S_kplus1 = 0.5*(S_kplus1+S_kplus1');
                NIS_vals(k+1) = innov_kplus1'*inv(S_kplus1)*innov_kplus1;
                
                
                estim_error_vals(k+1,:) = est_err_kplus1';

                
                break;   
                
            end
        end        
        
        % Log new estimate and est err covariance
        x_estim_vals(k+1, :) = x_kplus1';
        P_vals(:,:,k+1)      = P_kplus1;
        t_vals(k+1)          = t_kplus1;  
        
        
    end
    
    
    
    
    
%     % %         ---- Estimated state plots ----
%     figure()
%     suptitle('Extended Kalman Filter State Estimate');
%     set(findall(gcf,'type','text'),'FontSize',16)
% 
% 
%     subplot(4,1,1);
%     p1 = plot(t_vals, x_estim_vals(:,1), 'b', 'LineWidth',  0.1);
%     hold on;
%     plot(t_vals, x_Truth(:,1), 'b--', 'LineWidth',  0.5);
% %     plot(t_vals, x_nom_vals(:,1), 'b*', 'linewidth',  1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('x_{1}(t)','FontSize', 12)
%     grid on
% 
%     subplot(4,1,2)
%     p2 = plot(t_vals, x_estim_vals(:,2), 'r', 'LineWidth', 0.1);
%     hold on;
%     plot(t_vals, x_Truth(:,2), 'r--', 'LineWidth',  0.5);
% %     plot(t_vals, x_nom_vals(:,2), 'r*', 'linewidth',  1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('x_{2}(t)','FontSize', 12)
%     grid on
% 
%     subplot(4,1,3)
%     p3 = plot(t_vals, x_estim_vals(:,3), 'g', 'LineWidth', 0.1);
%     hold on;
%     plot(t_vals, x_Truth(:,3), 'g--', 'LineWidth',  0.1);
% %     plot(t_vals, x_nom_vals(:,3), 'g*', 'linewidth',  1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('x_{3}(t)','FontSize', 12)
%     grid on
% 
%     subplot(4,1,4)
%     p4 = plot(t_vals, x_estim_vals(:,4), 'k', 'LineWidth', 0.1);
%     hold on;
%     plot(t_vals, x_Truth(:,4), 'k--', 'LineWidth',  0.1);
% %     plot(t_vals, x_nom_vals(:,4), 'k*', 'linewidth',  1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('x_{4}(t)','FontSize', 12)
%     grid on
    
% %     % ---- Estimated Perturbation state plots ----
% %     figure()
% %     suptitle('Extended Kalman Filter State Estimate');
% %     set(findall(gcf,'type','text'),'FontSize',16)
% % 
% % 
% %     subplot(4,1,1);
% %     p1 = plot(t_vals, x_estim_vals(:,1), 'b-', 'linewidth',  1);
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$\hat{x_{1}}(t)$','Interpreter', 'latex', 'FontSize', 12)
% %     grid on
% % 
% %     subplot(4,1,2)
% %     p2 = plot(t_vals, x_estim_vals(:,2), 'r-', 'linewidth', 1);
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$\hat{x_{2}}(t)$','Interpreter', 'latex','FontSize', 12)
% %     grid on
% % 
% %     subplot(4,1,3)
% %     p3 = plot(t_vals, x_estim_vals(:,3), 'g-', 'linewidth', 1);
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$\hat{x_{3}}(t)$','Interpreter', 'latex','FontSize', 12)
% %     grid on
% % 
% %     subplot(4,1,4)
% %     p4 = plot(t_vals, x_estim_vals(:,4), 'k-', 'linewidth', 1);
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$\hat{x_{4}}(t)$','Interpreter', 'latex','FontSize', 12)
% %     grid on

end
