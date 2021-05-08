
% Run the LKF

function [x_pert_vals, x_estim_vals, P_vals, NEES_vals, NIS_vals] = ...
    RunLinearizedKalmanFilter(MEASUREMENTS_PROVIDED, y_Truth, ~, Q_true, R_true, x_Truth, ~, x_pert_0, P_0, num_steps)
%     RunLinearizedKalmanFilter(y_Truth, y_nom_vals, Q_true, R_true, x_Truth, x_nom_vals, x_pert_0, P_0, num_steps)

    global num_states delta_t xinit

%     % Check if the filter was provided measurement data, if not we'll have
%     % to simulate it ourselves
%     if any(class(y_Truth)) == 'cell'
%         MEASUREMENTS_PROVIDED = true;
%     else
%         MEASUREMENTS_PROVIDED = false;
%     end
    
    % Generate noisey nominal state information (TODO: should this include
    % noise?)
    tspan = 0:delta_t:num_steps*delta_t;            % Time span for ode45
    opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    u = zeros(2,1);
    w = zeros(2,1);
    [~, x_nom_vals] = ode45(@(t,y) nonLinearOde(t, y, u, w), tspan, xinit, opts);
    

    % Initialize arrays to store estimates and est err covariances
    x_estim_vals = zeros(num_steps+1, num_states);
    x_estim_vals(1, :) = (x_pert_0 + x_nom_vals(1,:));
    
    x_pert_vals = zeros(num_steps+1, num_states);
    x_pert_vals(1,:) = x_pert_0;
    
    t_vals = zeros(num_steps+1, 1);
    
    P_vals        = zeros(num_states, num_states, num_steps+1);
    P_vals(:,:,1) = P_0;
    
    NEES_vals    = zeros(num_steps+1,1);
    % If we have truth data available (currently, if MEASUREMENT_PROVIDED =
    % false, then truth data must be provided
    if not(isempty(x_Truth))
        est_err_0    = (x_Truth(1, :) - (x_nom_vals(1,:) + x_pert_0))'; % column vector    
        NEES_vals(1) = est_err_0'*inv(P_0)*est_err_0;
    end
    
    NIS_vals    = zeros(num_steps+1, 1);
    NIS_vals(1) = NaN;     % TODO: should this be NaN?                


    
% Tuning parameters
% HIGH Q -> low confidence in dynamics (trust the measurements)
% HIGH R -> low confidence in sensors (trust the prediction step)
%%Case 1: KF uses the correct/true Q and R values for DT system:     
%         Q_k      = 1*Q_true; 
        R_kplus1 = 1*R_true; 


%------------------Q1
%        Q_k      = [5*Q_true(1,1), 1e-11;
%                 1e-11, 5*Q_true(2,2)];               
%------------------
% % %        Q_k      = [Q_true(1,1), 1e-11;
% % %                 1e-11, Q_true(2,2)];
%-Q2
% % % %        Q_k      = [Q_true(1,1), 2e-11;
% % % %                 2e-11, Q_true(2,2)];
% -Q3
% % % %         Q_k      = [Q_true(1,1), 3e-11;
% % % %                 3e-11, Q_true(2,2)];
%-Q4
% % % % %        Q_k      = [2*Q_true(1,1), 3e-11;
% % % % %                 3e-11, 2*Q_true(2,2)];
%        Q_k      = [3*Q_true(1,1), 3e-11;
%                 3e-11, 3*Q_true(2,2)]; 
%        Q_k      = [0.75*Q_true(1,1), 3e-11;
%                 3e-11, 0.75*Q_true(2,2)];
%-Q5
       Q_k      = [2*Q_true(1,1), 5e-11;
                5e-11, 2*Q_true(2,2)];            



    for k =1:num_steps
    
      
        
               
        % Get values from t = tk
        x_nom_k      = x_nom_vals(k,:)';
%         Sw = chol(Q_true, 'lower');
%         w_k = Sw*(randn(2,1));
%         x_nom_k      = GetNominalState(k*delta_t) + [0 0; 1 0; 0 0; 0 1]*w_k;
%         x_nom_k      = GetNominalState(k*delta_t);
        dx_k         = x_pert_vals(k,:)';
        P_k          = P_vals(:,:,k);
        Omga_tilde_k = delta_t * [0 0; 1 0; 0 0; 0 1];  % Time invariant, Q = dT * Gamma(t)


        % --- Time update prediction step ---
        % Get linearized DT system parameters (evaluated at xnom, no noise or
        % inputs)
        unom = zeros(2,1);
        du_k = zeros(2,1);
        wprocess = zeros(2,1);
        [A_nom_eval, B_nom_eval] = LinearizedDynamicsOde(x_nom_k, unom, wprocess);

        % Compute the linearized DT system matrices
        F_tilde_k = eye(num_states, num_states) + delta_t * A_nom_eval; 
        G_tilde_k = delta_t * B_nom_eval;

        % Updated perturbation estimate and est error covariance
        dx_kplus1_min = F_tilde_k*dx_k + G_tilde_k*du_k;
        P_kplus1_min  = F_tilde_k*P_k*F_tilde_k' + Omga_tilde_k*Q_k*Omga_tilde_k';
        % ------------------------------------        

        % Get values from t = tk+1
        x_nom_kplus1    = x_nom_vals(k+1, :)';   % Nominal state at k+1
%         x_nom_kplus1 = GetNominalState((k+1)*delta_t);
        
        v_kplus1        = zeros(3,1);            % Process noise
        sim_time_kplus1 = (k+1) * delta_t;

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
                gsState = GetGroundStationState(sim_time_kplus1, gsIdx);
            end
            
        else
        % We have to simulate the measurement
        
            % Find the first usable measurement
            measurement_valid = false;
            gsIdx           = 1;                     % Initial ground station to get data from
            while not(measurement_valid)

                gsState = GetGroundStationState(sim_time_kplus1, gsIdx);

                % Get measurement from this ground station
                % Generate noise value to apply to measurement
                Sv = chol(R_true, 'lower');
                v = Sv*randn(3,1);
                y_tru_kplus1 = nonLinearMeasurementOde(sim_time_kplus1, x_Truth(k+1,:)', gsState, true) + v;

                % If the measurement is valid, use it and stop looking for
                % other measurements
                if not(isnan(y_tru_kplus1)) % && not(isnan(y_nom_kplus1)))

                    measurement_valid = true;

                % Otherwise, check the next ground station
                else
                    gsIdx = gsIdx + 1;
                end

                % If we have checked all the ground stations, move on
                if gsIdx > 12
                    break;
                end
            end
        
        end % end if MEASUREMENTS_PROVIDED
        
        
                
        % If we checked all ground stations and still have no data, just do
        % pure prediction
        if any(isnan(y_tru_kplus1)) || any(isempty(y_tru_kplus1))
            
            % No measurement update step, only pure prediction
            dx_kplus1 = dx_kplus1_min;
            P_kplus1  = P_kplus1_min;
            x_kplus1  = x_nom_kplus1 + dx_kplus1;
            
            % Don't evaluate NEES and NIS becuase our estimates are very
            % off
            NEES_vals(k+1) = NaN;
            NIS_vals(k+1) = NaN;
         
            
        else
            
            % --- Measurement update step --- 
            % Get the state and sensing matrix of the usable ground station
            H_tilde_kplus1 = LinearizedMeasurementOde(x_nom_kplus1, unom, v_kplus1, gsState);

            % Estimated measurement of nominal state
            y_nom_kplus1   = nonLinearMeasurementOde(sim_time_kplus1, x_nom_kplus1, gsState, false);

            
            % Compute Kalman gain
            K_kplus1 = P_kplus1_min*H_tilde_kplus1'*inv(H_tilde_kplus1*P_kplus1_min*H_tilde_kplus1' + R_kplus1);
           
            % Perturbation measurement
            dy_kplus1 = y_tru_kplus1 - y_nom_kplus1;

            % Updated perturbation estimate
            innov_kplus1 = dy_kplus1 - H_tilde_kplus1*dx_kplus1_min;

            % Wrap to pi
            innov_kplus1(end) = atan2(sin(innov_kplus1(end)), cos(innov_kplus1(end)));
            
            dx_kplus1 = dx_kplus1_min + K_kplus1*innov_kplus1;

            % Updated total state estimate
            x_kplus1 = x_nom_kplus1 + dx_kplus1;

            % Updated estimation error covariance
            P_kplus1 = (eye(num_states, num_states) - K_kplus1*H_tilde_kplus1)*P_kplus1_min;
            % ------------------------------------
            
            % If measurements were provided, we can't evaluate NEES and NIS
            if MEASUREMENTS_PROVIDED
                % skip
            else
                
                % Compute and Store NEES and NIS Statistics
                est_err_kplus1 = (x_Truth(k+1,:) - x_kplus1')';  % column vector
                NEES_vals(k+1) = est_err_kplus1'*inv(P_kplus1)*est_err_kplus1;

                S_kplus1 = H_tilde_kplus1*P_kplus1_min*H_tilde_kplus1' + R_kplus1;
                % DEBUGGING: using the same method as 1D robot example
                S_kplus1 = 0.5*(S_kplus1+S_kplus1');
                NIS_vals(k+1) = innov_kplus1'*inv(S_kplus1)*innov_kplus1;
            end
            
        end



        % Log new estimate and est err covariance
        x_pert_vals(k+1, :)  = dx_kplus1';
        x_estim_vals(k+1, :) = x_kplus1';
        P_vals(:,:,k+1)      = P_kplus1;
        t_vals(k+1)          = sim_time_kplus1;  
               
        
    end


% %     % ---- Estimated state plots ----
% %     figure()
% %     suptitle('Linearized Kalman Filter State Estimate');
% %     set(findall(gcf,'type','text'),'FontSize',18)
% % 
% % 
% %     subplot(4,1,1);
% %     p1 = plot(t_vals, x_estim_vals(:,1), 'b');
% %     hold on;
% %     plot(t_vals, x_Truth(:,1), 'b--');
% % %     plot(t_vals, x_nom_vals(:,1), 'b*', 'linewidth',  1);
% %     legend('Estimate', 'Truth');
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$x(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% %     grid on
% % 
% %     subplot(4,1,2)
% %     p2 = plot(t_vals, x_estim_vals(:,2), 'r');
% %     hold on;
% %     plot(t_vals, x_Truth(:,2), 'r--');
% % %     plot(t_vals, x_nom_vals(:,2), 'r*', 'linewidth',  1);
% %     legend('Estimate', 'Truth');
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$\dot{x}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% %     grid on
% % 
% %     subplot(4,1,3)
% %     p3 = plot(t_vals, x_estim_vals(:,3), 'g');
% %     hold on;
% %     plot(t_vals, x_Truth(:,3), 'g--');
% % %     plot(t_vals, x_nom_vals(:,3), 'g*', 'linewidth',  1);
% %     legend('Estimate', 'Truth');
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$y(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% %     grid on
% % 
% %     subplot(4,1,4)
% %     p4 = plot(t_vals, x_estim_vals(:,4), 'k');
% %     hold on;
% %     plot(t_vals, x_Truth(:,4), 'k--');
% % %     plot(t_vals, x_nom_vals(:,4), 'k*', 'linewidth',  1);
% %     legend('Estimate', 'Truth');
% %     xlabel('Time (s)','FontSize', 12)
% %     ylabel('$\dot{y}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% %     grid on
% %     set(findall(gcf,'type','line'),'linewidth',2)


% 
%     
%     % ---- Estimated Perturbation state plots ----
%     figure()
%     suptitle('Linearized Kalman Filter Perturbation State Estimate');
%     set(findall(gcf,'type','text'),'FontSize',16)
% 
% 
%     subplot(4,1,1);
%     p1 = plot(t_vals, x_pert_vals(:,1), 'b-', 'linewidth',  1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('\deltax_{1}(t)','FontSize', 12)
%     grid on
% 
%     subplot(4,1,2)
%     p2 = plot(t_vals, x_pert_vals(:,2), 'r-', 'linewidth', 1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('\deltax_{2}(t)','FontSize', 12)
%     grid on
% 
%     subplot(4,1,3)
%     p3 = plot(t_vals, x_pert_vals(:,3), 'g-', 'linewidth', 1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('\deltax_{3}(t)','FontSize', 12)
%     grid on
% 
%     subplot(4,1,4)
%     p4 = plot(t_vals, x_pert_vals(:,4), 'k-', 'linewidth', 1);
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('\deltax_{4}(t)','FontSize', 12)
%     grid on
    
    
    
end