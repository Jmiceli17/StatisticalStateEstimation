% EKF updated to simulate its own truth data


function [x_estim_vals, P_vals, NEES_vals, NIS_vals, estim_error_vals, x_Truth] = ...
    RunExtendedKalmanFilter_update(MEASUREMENTS_PROVIDED, Q_true, R_true, xinit, P_0, num_steps, y_Truth)

    global num_states delta_t r0 vel0 

    % Check if the filter was provided measurement data, if not we'll have
    % to simulate it ourselves
%     if (class(y_Truth) == 'cell')
%         MEASUREMENTS_PROVIDED = true;
%     else
%         MEASUREMENTS_PROVIDED = false;
%     end
    
    % ODE45 options and initialization
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    x_pert_0 = [0.0, 0.075, 0.0, -0.021]';          % Initial perturbed state
    x_nom    = [r0, 0.0, 0.0, vel0]';
    x_0_ode45 = x_nom + x_pert_0;
    
    x_0_ode45_truth = xinit;

    
    % Initialize arrays to store estimates and est err covariances
    x_estim_vals = zeros(num_steps+1, num_states);
    x_estim_vals(1, :) = xinit';
    
    t_vals = zeros(num_steps+1, 1);

    P_vals        = zeros(num_states, num_states, num_steps+1);
    P_vals(:,:,1) = P_0;
    
    NEES_vals    = zeros(num_steps+1,1);    
    NIS_vals    = zeros(num_steps+1, 1);
    NIS_vals(1) = NaN;     % TODO: should this be NaN?  

    estim_error_vals = zeros(num_steps+1, 4);
    x_Truth          = zeros(num_steps+1, 4);
    x_Truth(1,:)     = xinit';
    
    Q_k = Q_true;
    R_kplus1 = R_true;
%     R_kplus1 = 0.5*R_true;
%     R_kplus1 = 1.5*R_true;

% Q_k = 0.5*Q_true;

% Q_k = 1.85*Q_true;
Q_k = 5*Q_true;

% ----- Q1
%     Q_k = [Q_true(1,1), 1e-11;
%            1e-11, Q_true(2,2)];
%-------
%     Q_k = [Q_true(1,1), 5e-11;
%            5e-11, Q_true(2,2)];
%     Q_k = 0.1*Q_true;
%     Q_k = [0.9*Q_true(1,1), 1.5e-11;
%            1.5e-11, 0.9*Q_true(2,2)];


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
        [~, x_kplus1_min] = ode45(@(t,y) nonLinearOde(t, y, u, w), [0 delta_t], x_0_ode45, opts);  
        
        
        % Updated perturbation estimate and est error covariance
        x_kplus1_min  = x_kplus1_min(end,:)'; % last element from ode45 integration solution
        P_kplus1_min  = F_tilde_k*P_k*F_tilde_k' + Omga_tilde_k*Q_k*Omga_tilde_k';

        % Get values from t = tk+1
        t_kplus1 = delta_t*(k+1);
             
        
        if MEASUREMENTS_PROVIDED
            % Use the provided measurement (as a 4x1 array)
            y_tru_kplus1 = cell2mat(y_Truth(k+1));
                
            % If more than one measurement is provided, only take the first
            % one
            if size(y_tru_kplus1,2) > 1
                y_tru_kplus1 = y_tru_kplus1(:,1);
            end
            
            
            % Extract data from the provided measurement if it's not empty
            if not(any(isempty(y_tru_kplus1)))
                % The ground station index is provided as the last element of the measurement
                gsIdx = y_tru_kplus1(4);    

                % Trim the measurement
                y_tru_kplus1 = reshape(y_tru_kplus1(1:3), 3,1);

                % Get Ground station state
                gsState = GetGroundStationState(t_kplus1, gsIdx);
            end
            
        

            
        else
            % We have to simulate the measurement
            
            % Generate noisey truth data
            Sw = chol(Q_true, 'lower');
            wtilde_k = Sw*(randn(2,1));

            % Integrate the spacecraft dynamics with noise added as ZOH input
            [~, x_truth_kplus1] = ode45(@(t,y) nonLinearOde(t, y, u, wtilde_k), [0 delta_t], x_0_ode45_truth, opts);
            x_truth_kplus1  = x_truth_kplus1(end,:)';
            x_0_ode45_truth = x_truth_kplus1;
            
            
            % Get the first available measurement at t=tk+1
            measurement_valid = false;
            gsIdx = 1;
            while not(measurement_valid)
                
                
                % Get Ground station state
                gsState = GetGroundStationState(t_kplus1, gsIdx);

                % Generate noise value to apply to measurement
                Sv = chol(R_true, 'lower');
                v_k = Sv*randn(3,1);
                y_tru_kplus1 = nonLinearMeasurementOde(t_kplus1, x_truth_kplus1, gsState, true) + v_k;

                % If measurement estimate is valid, use it
                if any(isnan(y_tru_kplus1))
                    gsIdx = gsIdx + 1;
                else
                    measurement_valid = true;
                    break;
                end

                % If we have checked all the ground stations and still don't have a measurement, move on
                if gsIdx > 12
                    break;
                end
            end
            
        end
        
       

        % If measurement is not valid (or empty), skip
        if any(isnan(y_tru_kplus1)) || any(isempty(y_tru_kplus1))

            x_kplus1  = x_kplus1_min;
            P_kplus1  = P_kplus1_min;


            NEES_vals(k+1) = NaN;
            NIS_vals(k+1) = NaN;

        else

            
            
            % Compute gain
            H_tilde_kplus1 = LinearizedMeasurementOde(x_kplus1_min, u, v, gsState);
            K_kplus1       = P_kplus1_min*H_tilde_kplus1'*inv(H_tilde_kplus1*P_kplus1_min*H_tilde_kplus1' + R_kplus1);
            
            % Generate (noise-free) measurement estimate
            evaluate_visibility = false;
            y_kplus1_min = nonLinearMeasurementOde(t_kplus1, x_kplus1_min, gsState, evaluate_visibility );
            
            % Compute innovation vector
            innov_kplus1   = y_tru_kplus1 - y_kplus1_min;
            
            % Wrap to pi
            innov_kplus1(end) = atan2(sin(innov_kplus1(end)), cos(innov_kplus1(end)));

            % Updated total state estimate and estimate error covariance
            x_kplus1 = x_kplus1_min + K_kplus1*innov_kplus1;
            P_kplus1 = (eye(num_states, num_states) - K_kplus1*H_tilde_kplus1)*P_kplus1_min;

            % If measurements were provided, we can't evaluate NEES and NIS
            if MEASUREMENTS_PROVIDED
                % skip
            else
                % Compute and Store NEES and NIS Statistics
                est_err_kplus1 = x_truth_kplus1 - x_kplus1;  % column vector
                NEES_vals(k+1) = est_err_kplus1'*inv(P_kplus1)*est_err_kplus1;

                S_kplus1 = H_tilde_kplus1*P_kplus1_min*H_tilde_kplus1' + R_kplus1;

                % This guarantees that S is positive semi-definite (which it should be because it's a covariance matrix)
                S_kplus1 = 0.5*(S_kplus1+S_kplus1');
                NIS_vals(k+1) = innov_kplus1'*inv(S_kplus1)*innov_kplus1; 
           
                x_Truth(k+1,:)          = x_truth_kplus1';
                estim_error_vals(k+1,:) = est_err_kplus1';
                
            end
        end                                   
        
        % Initailize the integrator for the next loop
        x_0_ode45         = x_kplus1;
        
        % Log new estimate and est err covariance
        x_estim_vals(k+1, :)    = x_kplus1';
        P_vals(:,:,k+1)         = P_kplus1;
        t_vals(k+1)             = t_kplus1;

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
