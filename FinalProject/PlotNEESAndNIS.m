% PlotNEESAndNIS.m
%
% Description:
%   Plot NEES and NIS statistics from Monte-Carlo testing
%
% Inputs:
%   plot_title - string, title to be displayed on figure
%   NEES_TMT - double array of NEES values at each time step for each sim
%   alphaNEES - double bet. 0 and 1 (typically 0.05), significance interval
%               to be used in NEES testing
%   NIS_TMT - double array of NIS values at each time step for each sim
%   alphaNIS - double bet. 0 and 1 (typically 0.05), significance interval
%              to be used in NIS testing


function PlotNEESAndNIS(plot_title, NEES_TMT, alphaNEES, NIS_TMT, alphaNIS)

    global num_states

    num_simulations = size(NEES_TMT,2);
    
    epsBar_NEES = mean(NEES_TMT, 2); % returns the average of each row of NEES values (the average NEES at each time step)
    epsBar_NIS  = mean(NIS_TMT, 2);  % average NIS value at each time step

    % NEES confidence interval
    r1x_NEES = chi2inv(alphaNEES/2, num_simulations*num_states)./num_simulations;
    r2x_NEES = chi2inv(1-alphaNEES/2, num_simulations*num_states)./num_simulations;
    
    figure()
    subplot(2,1,1)
    plot(epsBar_NEES, 'ro', 'MarkerSize', 2, 'Linewidth', 2);
    grid on;
    hold on;
    plot(r1x_NEES*ones(size(epsBar_NEES)), 'r--', 'Linewidth', 2);
    plot(r2x_NEES*ones(size(epsBar_NEES)), 'r--', 'Linewidth', 2);
    NEES_title = 'NEES Estimation Results ' + plot_title;
    title(NEES_title,'FontSize',16)
    ylabel('$NEES, \bar{\epsilon}_x$','FontSize',14, 'Interpreter', 'latex')
    legend('NEES @ time k', 'r_1 bound', 'r_2 bound')
%     set(gca, 'YScale', 'log')
%     ylim([0,75]);

    
    
    % NIS confidence interval
    r1x_NIS = chi2inv(alphaNIS/2, num_simulations*3)./num_simulations;
    r2x_NIS = chi2inv(1-alphaNIS/2, num_simulations*3)./num_simulations;
    
    subplot(2,1,2)
    plot(epsBar_NIS, 'bo', 'MarkerSize', 2, 'Linewidth', 2);
    grid on;
    hold on;
    plot(r1x_NIS*ones(size(epsBar_NIS)), 'b--', 'Linewidth', 2);
    plot(r2x_NIS*ones(size(epsBar_NIS)), 'b--', 'Linewidth', 2);
    NIS_title = 'NIS Estimation Results ' + plot_title;
    title(NIS_title,'FontSize',16)
    xlabel('time step, k','FontSize',14)
    ylabel('$NIS, \bar{\epsilon}_y$','FontSize',14, 'Interpreter', 'latex')
    legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
%     set(gca, 'YScale', 'log')
%     ylim([0,20]);

    



end