% =========================================================================
% ASEN 5044 Statistical Estimation of Dynamical Systems
%
% Final Project
% Spacecraft Orbit Determination in 2D
%
% Joe Miceli 
% Jose Nido
%
% Description:
%   This final project involves simulating the linearized and and nonlinear
%   dynamics of a s/c in orbit. Ground stations provide measurements of the
%   s/c as range, range rate, and elevation angle. A linearized KF and
%   extended KF are implemented and evaluated using Monte-Carlo truth model
%   testing. Note that some of the plotting sections may be currently
%   commented out.
% =========================================================================
clc; clear; close all;

% Constants and globals used through out the script
global mu r0 re vel0 xinit num_states num_inputs delta_t orb_period n xnom we

mu         = 398600.0;   % km^3/s^2
r0         = 6678.0;     % Nominal orbit radius, km
re         = 6378.0;     % Radius of earth, km
we         = 2*pi/86400; % rad/s
vel0       = r0*sqrt(mu/(r0^3));
xinit      = [r0, 0.0, 0.0, vel0]';
num_states = 4;
num_inputs = 2;
delta_t    = 10; % sec
orb_period = 2*pi*sqrt(r0^3/mu);    % seconds
n = 2*pi/orb_period;    % mean motion, rad/sec
% Symbolic nominal state
syms t
xnom = [r0*cos(n*t); -vel0*sin(n*t); r0*sin(n*t); vel0*cos(n*t)];



%%%%%%%%%%%%%%%%%% PART 1: DETERMINISITC SYSTEM ANALYSIS %%%%%%%%%%%%%%%%%%
% Define the state as X = [xpos, xvel, ypos, yvel]

%%% [Problem 1] 
% Compute CT Jacobians to obtain linearized parameters


%%% [Problem 2] 
% Discretize the linearized system (assuming no process 
% noise for now)
% Becuase the system is time varying, the discretized matrix changes each
% time step. Therefore, an observability and controllability anlaysis for
% the system cannot be conducted.


%%% [Problem 3]
% Simulate the linearized model and compare it to simulated non-linear
% model (using numerical integration scheme)
x_pert_0 = [0.0, 0.075, 0.0, -0.021]';          % Initial perturbed state
u_pert_0 = zeros(2,1);                          % Initial perturbed input
num_steps = 1400;                               
tspan = 0:delta_t:num_steps*delta_t;            % Time span for ode45
num_GrndStations = 12;                          % # ground stations
x0 = (xinit + x_pert_0)';                       % Propagate x0 + perturb_x0

% Run the simulation using the linearized dynamics (no noise)
[t_vals_LZ, x_LZ, xnom_evals_LZ, y_nom_vals_LZ] = SimulateLinearizedSystem(x_pert_0, u_pert_0, num_steps, num_GrndStations);

% Run the simulation using the nonlinear dynamics (no noise)
[t_vals_NL, x_NL, y_NL] = SimulateNonlinearSystem(tspan, x0, num_GrndStations);

% % % ---- Linearized and Nonlinear State Plot ----
% % figure()
% % suptitle('Simulated State, Linearized vs Nonlinear Integration' );
% % set(findall(gcf,'type','text'),'FontSize',18)
% % 
% % 
% % subplot(4,1,1);
% % p1 = plot(t_vals_LZ, x_LZ(:,1), 'b-', 'linewidth',  2);
% % hold on;
% % plot(t_vals_NL, x_NL(:,1), 'b--', 'linewidth',  2);
% % legend('Linearized', 'Non-linear')
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$x(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% % grid on
% % 
% % subplot(4,1,2)
% % p2 = plot(t_vals_LZ, x_LZ(:,2), 'r-', 'linewidth', 2);
% % hold on;
% % plot(t_vals_NL, x_NL(:,2), 'r--', 'linewidth',  2);
% % legend('Linearized', 'Non-linear')
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\dot{x}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % grid on
% % 
% % subplot(4,1,3)
% % p3 = plot(t_vals_LZ, x_LZ(:,3), 'g-', 'linewidth', 2);
% % hold on;
% % plot(t_vals_NL, x_NL(:,3), 'g--', 'linewidth',  2);
% % legend('Linearized', 'Non-linear')
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$y(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% % grid on
% % 
% % subplot(4,1,4)
% % p4 = plot(t_vals_LZ, x_LZ(:,4), 'k-', 'linewidth', 2);
% % hold on;
% % plot(t_vals_NL, x_NL(:,4), 'k--', 'linewidth',  2);
% % legend('Linearized', 'Non-linear')
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\dot{y}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % grid on
% % set(findall(gcf,'type','line'),'linewidth',2)


%%%%%%%%%%%%%%%%% PART 2: STOCHASTIC NONLINEAR FILTERING %%%%%%%%%%%%%%%%%
load('orbitdeterm_finalproj_KFdata.mat')

num_simulations = 10;


% P_0 = 0.01*eye(num_states,num_states);
% P_0(2,2) = 0.0001;
% P_0(4,4) = 0.0001;
P_0 = 0.1*eye(num_states,num_states);
P_0(2,2) = 0.001;
P_0(4,4) = 0.001;

x_pert_0_test = zeros(4,1);

useProcNoise = 1;



% --- Run Truth model testing ---

% Arrays to hold truth data outputs
x_truth_TMT          = zeros(num_steps+1, num_states, num_simulations);
y_truth_TMT          = zeros(num_steps+1, 3, num_GrndStations, num_simulations);
visible_gs_array_TMT = zeros(num_steps+1, num_GrndStations, num_simulations);

% Arrays for to hold LKF outputs
x_pert_estim_TMT = zeros(num_steps+1, num_states, num_simulations);
x_estim_TMT      = zeros(num_steps+1, num_states, num_simulations);
P_TMT            = zeros(num_states, num_states, num_steps+1, num_simulations);
NEES_TMT         = zeros(num_steps + 1, num_simulations);
NIS_TMT          = zeros(num_steps + 1, num_simulations);
estim_error_TMT      = zeros(num_steps + 1, num_states, num_simulations);
sigma1_TMT       = zeros(num_steps + 1, num_simulations);
sigma2_TMT       = zeros(num_steps + 1, num_simulations);
sigma3_TMT       = zeros(num_steps + 1, num_simulations);
sigma4_TMT       = zeros(num_steps + 1, num_simulations);


% Arrays for to hold EKF outputs
x_estim_TMT_EKF      = zeros(num_steps+1, num_states, num_simulations);
P_TMT_EKF            = zeros(num_states, num_states, num_steps+1, num_simulations);
NEES_TMT_EKF         = zeros(num_steps + 1, num_simulations);
NIS_TMT_EKF          = zeros(num_steps + 1, num_simulations);
estim_error_TMT_EKF  = zeros(num_steps + 1, num_states, num_simulations);
sigma1_TMT_EKF       = zeros(num_steps + 1, num_simulations);
sigma2_TMT_EKF       = zeros(num_steps + 1, num_simulations);
sigma3_TMT_EKF       = zeros(num_steps + 1, num_simulations);
sigma4_TMT_EKF       = zeros(num_steps + 1, num_simulations);



for ss = 1:num_simulations
    
    % Generate truth data 
    [t_truth_vals, x_truth_vals, y_truth_vals, y_nom_truth_vals, visible_gs_array_truth] =...
    GenerateTruthData(useProcNoise, num_steps, num_GrndStations, xinit, Qtrue, Rtrue);
    
    % Log truth data
    x_truth_TMT(:,:,ss)          = x_truth_vals;
    y_truth_TMT(:,:,:,ss)        = y_truth_vals;
    visible_gs_array_TMT(:,:,ss) = visible_gs_array_truth;

    % Run the linearized kalman filter 
    external_measurements = false; % we're using measurements we generated
    [x_pert_estim_vals, x_estim_vals, P_vals, NEES_vals, NIS_vals] =...
        RunLinearizedKalmanFilter(external_measurements, y_truth_vals, y_nom_truth_vals, Qtrue,...
        Rtrue, x_truth_vals, xnom_evals_LZ, x_pert_0_test', P_0, num_steps);
    
    % Log KF data
    x_pert_estim_TMT(:,:,ss) = x_pert_estim_vals;
    x_estim_TMT(:,:,ss)      = x_estim_vals;
    estim_error_TMT(:,:,ss)  = x_truth_vals - x_estim_vals;
    P_TMT(:,:,:,ss)          = P_vals;
    sigma1_TMT(:,ss)         = sqrt(P_vals(1,1,:));
    sigma2_TMT(:,ss)         = sqrt(P_vals(2,2,:));
    sigma3_TMT(:,ss)         = sqrt(P_vals(3,3,:));
    sigma4_TMT(:,ss)         = sqrt(P_vals(4,4,:));
    
    NEES_TMT(:,ss)           = NEES_vals; 
    NIS_TMT(:,ss)            = NIS_vals;
    
    
    
    
%     % Run the extended kalman filter 
%    [x_estim_vals_EKF, P_vals_EKF, NEES_vals_EKF, NIS_vals_EKF, estim_error_vals] = ...
%     RunExtendedKalmanFilter(y_truth_vals, Qtrue, Rtrue, x_truth_vals, xinit, P_0, num_steps);

%     % Run the extended KF v2.0
     [x_estim_vals_EKF, P_vals_EKF, NEES_vals_EKF, NIS_vals_EKF, estim_error_vals, x_truth_update] =...
    RunExtendedKalmanFilter_update(external_measurements, Qtrue, Rtrue, xinit, P_0, num_steps, []);


    % Log KF data
    x_estim_TMT_EKF(:,:,ss)      = x_estim_vals_EKF;
    estim_error_TMT_EKF(:,:,ss)  = estim_error_vals;
    P_TMT_EKF(:,:,:,ss)          = P_vals_EKF;
    sigma1_TMT_EKF(:,ss)         = sqrt(P_vals_EKF(1,1,:));
    sigma2_TMT_EKF(:,ss)         = sqrt(P_vals_EKF(2,2,:));
    sigma3_TMT_EKF(:,ss)         = sqrt(P_vals_EKF(3,3,:));
    sigma4_TMT_EKF(:,ss)         = sqrt(P_vals_EKF(4,4,:));
    
    NEES_TMT_EKF(:,ss)           = NEES_vals_EKF; 
    NIS_TMT_EKF(:,ss)            = NIS_vals_EKF;
    
end 
% % % % % 
% % % % % % Significance levels
% % % % % alphaNEES = 0.01;
% % % % % alphaNIS = 0.01;
% % % % % 
% % % % % LKF_title = "Linearized Kalman Filter";
% % % % % PlotNEESAndNIS(LKF_title, NEES_TMT, alphaNEES, NIS_TMT, alphaNIS)
% % % % % 
% % % % % EKF_title = "Extended Kalman Filter";
% % % % % PlotNEESAndNIS(EKF_title, NEES_TMT_EKF, alphaNEES, NIS_TMT_EKF, alphaNIS)
% % % % % 
% % % % % 
% % % % % 
% % % % % plot_title_LKF = 'Typical Estimation Error and +2\sigma Bounds, Linearized Kalman Filter';
% % % % % PlotErrorAndSigma(plot_title_LKF, t_truth_vals, sigma1_TMT(:,end), sigma2_TMT(:,end), sigma3_TMT(:,end), sigma4_TMT(:,end), estim_error_TMT(:,:,end))
% % % % % 
% % % % % 
% % % % % plot_title_EKF = 'Typical Estimation Error and +2\sigma Bounds, Extended Kalman Filter';
% % % % % PlotErrorAndSigma(plot_title_EKF, t_truth_vals, sigma1_TMT_EKF(:,end), sigma2_TMT_EKF(:,end), sigma3_TMT_EKF(:,end), sigma4_TMT_EKF(:,end), estim_error_TMT_EKF(:,:,end))
% % % % % 
% % % % % 
% % % % % % plot_title_EKF = 'Typical Estimation Error and +2\sigma Bounds, Extended Kalman Filter';
% % % % % % PlotErrorAndSigma(plot_title_EKF, t_truth_vals, sigma1_TMT_EKF(:,end), sigma2_TMT_EKF(:,end), sigma3_TMT_EKF(:,end), sigma4_TMT_EKF(:,end), estim_error_vals)



% % % --- Plot 2sigma ---
% % figure()
% % suptitle(plot_title);
% % set(findall(gcf,'type','text'),'FontSize',16)
% % 
% % 
% % 
% % subplot(4,1,1);
% % p1 = plot(t_truth_vals, 2*abs(sigma1_TMT(:,:,end)), 'k--', 'linewidth',  2);
% % hold on;
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma_{11} [km^{2}]$', 'Interpreter', 'latex','FontSize', 12)
% % grid on
% % 
% % subplot(4,1,2)
% % p2 = plot(t_truth_vals, 2*abs(sigma2_TMT(:,:,end)), 'k--', 'linewidth', 2);
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma{22} [km^{2}/s^{2}]$', 'Interpreter', 'latex', 'FontSize', 12)
% % grid on
% % 
% % subplot(4,1,3)
% % p3 = plot(t_truth_vals, 2*abs(sigma3_TMT(:,:,end)), 'k--', 'linewidth', 2);
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma_{33} [km^{2}]$', 'Interpreter', 'latex', 'FontSize', 12)
% % grid on
% % 
% % subplot(4,1,4)
% % p4 = plot(t_truth_vals, 2*abs(sigma4_TMT(:,:,end)), 'k--', 'linewidth', 2);
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma_{44} [km^{2}/s^{2}]$', 'Interpreter', 'latex','FontSize', 12)
% % grid on


% % % --- Plot 2sigma EKF ---
% % figure()
% % suptitle(plot_title_EKF);
% % set(findall(gcf,'type','text'),'FontSize',16)
% % 
% % 
% % 
% % subplot(4,1,1);
% % p1 = plot(t_truth_vals, 2*abs(sigma1_TMT_EKF(:,end)), 'k--', 'linewidth',  2);
% % hold on;
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma_{11} [km^{2}]$', 'Interpreter', 'latex','FontSize', 12)
% % grid on
% % 
% % subplot(4,1,2)
% % p2 = plot(t_truth_vals, 2*abs(sigma2_TMT_EKF(:,end)), 'k--', 'linewidth', 2);
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma{22} [km^{2}/s^{2}]$', 'Interpreter', 'latex', 'FontSize', 12)
% % grid on
% % 
% % subplot(4,1,3)
% % p3 = plot(t_truth_vals, 2*abs(sigma3_TMT_EKF(:,end)), 'k--', 'linewidth', 2);
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma_{33} [km^{2}]$', 'Interpreter', 'latex', 'FontSize', 12)
% % grid on
% % 
% % subplot(4,1,4)
% % p4 = plot(t_truth_vals, 2*abs(sigma4_TMT_EKF(:,end)), 'k--', 'linewidth', 2);
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\sigma_{44} [km^{2}/s^{2}]$', 'Interpreter', 'latex','FontSize', 12)
% % grid on


% % % ---- Truth value plots ----
% % figure()
% % suptitle('Typical Simulation State Truth and Linearized KF Estimate');
% % set(findall(gcf,'type','text'),'FontSize',16)
% % 
% % 
% % subplot(4,1,1);
% % p1 = plot(t_truth_vals, x_truth_vals(:,1), 'b--', 'linewidth',  2);
% % hold on;
% % plot(t_truth_vals, x_estim_vals(:,1), 'b-', 'linewidth',  2);
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('x_{1}(t) [km]','FontSize', 12)
% % grid on
% % 
% % subplot(4,1,2)
% % p2 = plot(t_truth_vals, x_truth_vals(:,2), 'r--', 'linewidth', 2);
% % hold on;
% % plot(t_truth_vals, x_estim_vals(:,2), 'r-', 'linewidth',  2);
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('x_{2}(t) [km/s]','FontSize', 12)
% % grid on
% % 
% % subplot(4,1,3)
% % p3 = plot(t_truth_vals, x_truth_vals(:,3), 'g--', 'linewidth', 2);
% % hold on;
% % plot(t_truth_vals, x_estim_vals(:,3), 'g-', 'linewidth',  2);
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('x_{3}(t) [km]','FontSize', 12)
% % grid on
% % 
% % subplot(4,1,4)
% % p4 = plot(t_truth_vals, x_truth_vals(:,4), 'k--', 'linewidth', 2);
% % hold on;
% % plot(t_truth_vals, x_estim_vals(:,4), 'k-', 'linewidth',  2);
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('x_{4}(t) [km/s]','FontSize', 12)
% % grid on
% % 




% % % ---- EKF Truth value plots ----
% % figure()
% % suptitle('Typical Simulation State Truth and EKF Estimate OLD!!!');
% % set(findall(gcf,'type','text'),'FontSize',16)
% % 
% % 
% % subplot(4,1,1);
% % p1 = plot(t_truth_vals, x_truth_vals(:,1), 'b--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,1), 'b-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\hat{x_{1}}(t) [km]$','Interpreter', 'latex', 'FontSize', 12)
% % grid on
% % 
% % subplot(4,1,2)
% % p2 = plot(t_truth_vals, x_truth_vals(:,2), 'r--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,2), 'r-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\hat{x_{2}}(t) [km/s]$','Interpreter', 'latex','FontSize', 12)
% % grid on
% % 
% % subplot(4,1,3)
% % p3 = plot(t_truth_vals, x_truth_vals(:,3), 'g--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,3), 'g-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\hat{x_{3}}(t) [km]$','Interpreter', 'latex','FontSize', 12)
% % grid on
% % 
% % subplot(4,1,4)
% % p4 = plot(t_truth_vals, x_truth_vals(:,4), 'k--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,4), 'k-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\hat{x_{4}}(t) [km/s]$','Interpreter', 'latex','FontSize', 12)
% % grid on
% % 
% % 
% % % ---- EKF Truth value plots ----
% % figure()
% % suptitle('Typical Simulation State Truth and EKF Estimate');
% % set(findall(gcf,'type','text'),'FontSize',18)
% % 
% % 
% % subplot(4,1,1);
% % p1 = plot(t_truth_vals, x_truth_update(:,1), 'b--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,1), 'b-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$x(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% % grid on
% % 
% % subplot(4,1,2)
% % p2 = plot(t_truth_vals, x_truth_update(:,2), 'r--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,2), 'r-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\dot{x}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % grid on
% % 
% % subplot(4,1,3)
% % p3 = plot(t_truth_vals, x_truth_update(:,3), 'g--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,3), 'g-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$y(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% % grid on
% % 
% % subplot(4,1,4)
% % p4 = plot(t_truth_vals, x_truth_update(:,4), 'k--');
% % hold on;
% % plot(t_truth_vals, x_estim_vals_EKF(:,4), 'k-');
% % legend('Truth', 'Estimate');
% % xlabel('Time (s)','FontSize', 12)
% % ylabel('$\dot{y}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% % grid on
% % set(findall(gcf,'type','line'),'linewidth',2)




% % ---- Measurement State Plot ----
% figure()
% suptitle('Full Noisy Nonlinear Model Data Simulation');
% set(findall(gcf,'type','text'),'FontSize',16)
% 
% 
% for gsIdx = 1:num_GrndStations
%     
%     % Generate random color scheme to use for this ground station's data
%     rand_color = [rand, rand, rand];
%     
%     subplot(4,1,1);
%     hold on;
%     plot(t_truth_vals, y_truth_vals(:,1,gsIdx), '-', 'Color', rand_color);
%     hold off;
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('\rho(t) [km]','FontSize', 12)
%     grid on
%     
%     subplot(4,1,2)
%     hold on;
%     plot(t_truth_vals, y_truth_vals(:,2,gsIdx), '-', 'Color', rand_color);
%     hold off;
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('$\dot{\rho}(t) [km/s]$','FontSize', 12, 'Interpreter', 'latex')
%     grid on
%     
%     subplot(4,1,3)
%     hold on;
%     plot(t_truth_vals, y_truth_vals(:,3,gsIdx), '-', 'Color', rand_color);
%     hold off;
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('\phi(t) [rad]','FontSize', 12)
%     grid on
%     
%     subplot(4,1,4)
%     hold on;
%     plot(t_truth_vals, visible_gs_array_truth(:,gsIdx), '*', 'Color', rand_color);
%     hold off;
%     xlabel('Time (s)','FontSize', 12)
%     ylabel('Visible Ground Station ID','FontSize', 12)
%     grid on
%     
% end



%%% [Problem 6]
% Estimate the state of the satellite using the LKF and EKF with the 
% provided measurement data

x_init_P6 = xinit;
P_0_P6 = diag([0.1, 0.001, 0.1, 0.001]);


% Run the linearized kalman filter 
meas_provided = true;

[~, x_estim_vals_LKF_P6, P_vals_LKF_P6, ~, ~] =...
    RunLinearizedKalmanFilter(meas_provided, ydata, [], Qtrue,...
    Rtrue, [], xnom_evals_LZ, x_pert_0_test', P_0, num_steps);

% variance values from LKF
sigma1_LKF_P6 = reshape(sqrt(P_vals_LKF_P6(1,1,:)), size(P_vals_LKF_P6,3), 1);
sigma2_LKF_P6 = reshape(sqrt(P_vals_LKF_P6(2,2,:)), size(P_vals_LKF_P6,3), 1);
sigma3_LKF_P6 = reshape(sqrt(P_vals_LKF_P6(3,3,:)), size(P_vals_LKF_P6,3), 1);
sigma4_LKF_P6 = reshape(sqrt(P_vals_LKF_P6(4,4,:)), size(P_vals_LKF_P6,3), 1);


% Run the extended KF v2.0
 [x_estim_vals_EKF_P6, P_vals_EKF_P6, ~, ~, ~, ~] =...
RunExtendedKalmanFilter_update(meas_provided, Qtrue, Rtrue, xinit, P_0_P6, num_steps, ydata);

% variance values from EKF
sigma1_EKF_P6 = reshape(sqrt(P_vals_EKF_P6(1,1,:)), size(P_vals_EKF_P6,3), 1);
sigma2_EKF_P6 = reshape(sqrt(P_vals_EKF_P6(2,2,:)), size(P_vals_EKF_P6,3), 1);
sigma3_EKF_P6 = reshape(sqrt(P_vals_EKF_P6(3,3,:)), size(P_vals_EKF_P6,3), 1);
sigma4_EKF_P6 = reshape(sqrt(P_vals_EKF_P6(4,4,:)), size(P_vals_EKF_P6,3), 1);


% % ---- EKF & LKF Truth value plots ----
% figure()
% suptitle('LKF and EKF State Estimate Using Provided Observations');
% set(findall(gcf,'type','text'),'FontSize',18)
% 
% 
% subplot(4,1,1);
% plot(t_vals_NL, x_estim_vals_EKF_P6(:,1), 'b-');
% hold on;
% plot(t_vals_NL, x_estim_vals_LKF_P6(:,1), 'b--');
% legend({'EKF Estimate','LKF Estimate'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$x(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% grid on
% 
% subplot(4,1,2)
% plot(t_vals_NL, x_estim_vals_EKF_P6(:,2), 'r-');
% hold on;
% plot(t_vals_NL, x_estim_vals_LKF_P6(:,2), 'r--');
% legend({'EKF Estimate','LKF Estimate'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$\dot{x}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% grid on
% 
% subplot(4,1,3)
% plot(t_vals_NL, x_estim_vals_EKF_P6(:,3), 'g-');
% hold on;
% plot(t_vals_NL, x_estim_vals_LKF_P6(:,3), 'g--');
% legend({'EKF Estimate','LKF Estimate'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$y(t) [km]$','FontSize', 14, 'Interpreter', 'latex')
% grid on
% 
% subplot(4,1,4)
% plot(t_vals_NL, x_estim_vals_EKF_P6(:,4), 'k-');
% hold on;
% plot(t_vals_NL, x_estim_vals_LKF_P6(:,4), 'k--');
% legend({'EKF Estimate','LKF Estimate'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$\dot{y}(t) [km/s]$','FontSize', 14,'Interpreter', 'latex')
% grid on
% set(findall(gcf,'type','line'),'linewidth',2)
% 
% 
% 
% 
% % --- Plot 2sigma EKF & LKF ---
% figure()
% suptitle('LKF and EKF 2\sigma Estimation Error Variance');
% set(findall(gcf,'type','text'),'FontSize',18)
% 
% subplot(4,1,1);
% plot(t_vals_NL, 2*abs(sigma1_LKF_P6(:,end)), 'k-', 'linewidth',  1.5);
% hold on;
% plot(t_vals_NL, 2*abs(sigma1_EKF_P6(:,end)), 'r-', 'linewidth',  1.5);
% legend({'LKF','EKF'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$\sigma_{11} [km]$', 'Interpreter', 'latex','FontSize', 14)
% grid on
% 
% subplot(4,1,2)
% plot(t_vals_NL, 2*abs(sigma2_LKF_P6(:,end)), 'k-', 'linewidth',  1.5);
% hold on;
% plot(t_vals_NL, 2*abs(sigma2_EKF_P6(:,end)), 'r-', 'linewidth', 1.5);
% legend({'LKF','EKF'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$\sigma_{22} [km/s]$', 'Interpreter', 'latex', 'FontSize', 14)
% grid on
% 
% subplot(4,1,3)
% plot(t_vals_NL, 2*abs(sigma3_LKF_P6(:,end)), 'k-', 'linewidth',  1.5);
% hold on;
% plot(t_vals_NL, 2*abs(sigma3_EKF_P6(:,end)), 'r-', 'linewidth', 1.5);
% legend({'LKF','EKF'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$\sigma_{33} [km]$', 'Interpreter', 'latex', 'FontSize', 14)
% grid on
% 
% subplot(4,1,4)
% plot(t_vals_NL, 2*abs(sigma4_LKF_P6(:,end)), 'k-', 'linewidth',  1.5);
% hold on;
% p4 = plot(t_vals_NL, 2*abs(sigma4_EKF_P6(:,end)), 'r-', 'linewidth', 1.5);
% legend({'LKF','EKF'}, 'Fontsize', 14);
% xlabel('Time (s)','FontSize', 12)
% ylabel('$\sigma_{44} [km/s]$', 'Interpreter', 'latex','FontSize', 14)
% grid on

