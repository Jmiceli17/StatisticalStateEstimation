% PlotErrorAndSigma.m
%
% Description:
%   Plot estimation error and 2sigma bounds from a filter


function PlotErrorAndSigma(plot_title, t_truth_vals, sigma1, sigma2, sigma3, sigma4, estim_error)

% ---- Sigma plots of last simulation ----
figure()
suptitle(plot_title);
set(findall(gcf,'type','text'),'FontSize',18)



subplot(4,1,1);
p1 = plot(t_truth_vals, 2*sigma1, 'k--', 'linewidth',  2);
hold on;
plot(t_truth_vals, estim_error(:,1), 'b-', 'linewidth',  2);
plot(t_truth_vals, -2*sigma1, 'k--', 'linewidth',  2);
xlabel('Time (s)','FontSize', 12)
ylabel('$e_{x} [km]$', 'Interpreter', 'latex','FontSize', 14)
grid on

subplot(4,1,2)
p2 = plot(t_truth_vals, 2*sigma2, 'k--', 'linewidth', 2);
hold on;
plot(t_truth_vals, estim_error(:,2), 'r-', 'linewidth',  2);
plot(t_truth_vals, -2*sigma2, 'k--', 'linewidth',  2);
xlabel('Time (s)','FontSize', 12)
ylabel('$e_{\dot{x}} [km/s]$', 'Interpreter', 'latex', 'FontSize', 14)
grid on

subplot(4,1,3)
p3 = plot(t_truth_vals, 2*sigma3, 'k--', 'linewidth', 2);
hold on;
plot(t_truth_vals, estim_error(:,3), 'g-', 'linewidth',  2);
plot(t_truth_vals, -2*sigma3, 'k--', 'linewidth',  2);
xlabel('Time (s)','FontSize', 12)
ylabel('$e_{y} [km]$', 'Interpreter', 'latex', 'FontSize', 14)
grid on

subplot(4,1,4)
p4 = plot(t_truth_vals, 2*sigma4, 'k--', 'linewidth', 2);
hold on;
plot(t_truth_vals, estim_error(:,4), 'k-', 'linewidth',  2);
plot(t_truth_vals, -2*sigma4, 'k--', 'linewidth',  2);
xlabel('Time (s)','FontSize', 12)
ylabel('$e_{\dot{y}} [km/s]$', 'Interpreter', 'latex', 'FontSize', 14)
grid on
set(findall(gcf,'type','line'),'linewidth',2)






end
