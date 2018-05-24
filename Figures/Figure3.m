%% Figure 3
% This script reproduces Figure 3 in:
%
% Wallot, S., Roepstorff, A., & Mønster, D. (2016). Multidimensional
% Recurrence Quantification Analysis (MdRQA) for the analysis of
% multidimensional time-series: A software implementation in MATLAB and its
% application to group-level data in joint action. Frontiers in Psychology,
% 7, 1835. http://dx.doi.org/10.3389/fpsyg.2016.01835
%
% Depends on:
% -Signal Processing Toolbox
% -Cross Recurrence Plot Toolbox
% -MdRQA script (MDRQA.m)
% -lorenz_model.m
% -set_fonts.m
%
% Description:
% RQA and MDRQA analysis of the Lorenz model.
% This script applies Recurrence Quantification Analysis (RQA) and
% Multi-Dimensional Recurrence Quantification Analysis (MDRQA) to the
% Lorenz model.
%
% By:
% Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany
% Dan Mønster, Aarhus University, Aarhus, Denmark
%
% July 2016.
%

fontsize = 14; % Font size to use in plots.

% ===============================
% Parameters for the Lorenz model
% ===============================
x0 = 10;
y0 = 10;
z0 = 10;
sigma = 10;
rho = 28;
beta = 8/3.0;
tmax = 20;

% ===============================
% Perform numerical integration
% ===============================
[t, x, y, z] = lorenz_model([x0, y0 z0], sigma, rho, beta, tmax);
xx = zscore(x.Data);
yy = zscore(y.Data);
zz = zscore(z.Data);

%% Plot the time series
figure('units','normalized','position',[.1 .1 .7 .9])
subplot(3,3,1)
plot(x,'k','linewidth',1);
xlabel('$t$','interpreter','latex');
ylabel('$x(t)$','interpreter','latex');
title('');
set_fonts(gcf, gca, fontsize);

subplot(3,3,2)
plot(y,'k','linewidth',1);
xlabel('$t$','interpreter','latex');
ylabel('$y(t)$','interpreter','latex');
title('');
set_fonts(gcf, gca, fontsize);

subplot(3,3,3)
plot(z,'k','linewidth',1);
xlabel('$t$','interpreter','latex');
ylabel('$z(t)$','interpreter','latex');
title('');
set_fonts(gcf, gca, fontsize);


%% Plot the attractor
subplot(3,4,4+4)
plot3(xx, yy, zz,'k');
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
zlabel('$z$','interpreter','latex')
xlim([-2.5, 2.5]);
ylim([-2.5, 2.5]);
zlim([-2, 2.5]);
view([45 30]);
set_fonts(gcf, gca, fontsize);

%% Resample the time series
% To get a uniform sampling that is the same for all three time series
Nsamples = numel(t);
sampling_interval = tmax/Nsamples;
[tx,  ty] = synchronize(x,y,'Uniform','Interval',sampling_interval);
[t2x, tz] = synchronize(x,z,'Uniform','Interval',sampling_interval);
if ~all(tx.Time == t2x.Time)
    disp('Problem synchronizing')
end
xx = zscore(tx.Data);
yy = zscore(ty.Data);
zz = zscore(tz.Data);

%% Set phase-space embedding parameters
T = 3;
tau = 4;
% RQA radius
radius = 0.1;

%% Reconstruct the attractor using x
subplot(3,4,4+1)
M = psembed(xx, T, tau);
plot3(M(:,1), M(:,2), M(:,3), 'k');
xlabel('$\tilde{V}_1(x)$','interpreter','latex')
ylabel('$\tilde{V}_2(x)$','interpreter','latex')
zlabel('$\tilde{V}_3(x)$','interpreter','latex')
xlim([-2.2, 2.2]);
ylim([-2.5, 2.5]);
zlim([-2.5, 2.5]);
view([45 30]);
set_fonts(gcf, gca, fontsize);

%% RP using x reconstruction
subplot(3,4,4+5)

rp = crp(xx, xx, T, tau,radius,...
    'euclidean','nonormalize','silent');
rqa_x = crqa(xx, xx, T, tau,radius,...
    'euclidean','nonormalize','silent');
spy(rp,'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set_fonts(gcf, gca, fontsize);


%% Reconstruct the attractor using y
subplot(3,4,4+2)
M = psembed(yy, T, tau);
plot3(M(:,1), M(:,2), M(:,3), 'k');
xlabel('$\tilde{V}_1(y)$','interpreter','latex')
ylabel('$\tilde{V}_2(y)$','interpreter','latex')
zlabel('$\tilde{V}_3(y)$','interpreter','latex')
xlim([-2.5, 2.5]);
ylim([-2.5, 2.5]);
zlim([-2.5, 2.5]);
view([45 30]);
set_fonts(gcf, gca, fontsize);

%% RP using y reconstruction
subplot(3,4,4+6)

rp = crp(yy, yy, T, tau,radius,...
    'euclidean','nonormalize','silent');
rqa_y = crqa(yy, yy, T, tau,radius,...
    'euclidean','nonormalize','silent');
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set_fonts(gcf, gca, fontsize);


%% Reconstruct the attractor using z
subplot(3,4,4+3)
M = psembed(zz, T, tau);
plot3(M(:,1), M(:,2), M(:,3), 'k');
xlabel('$\tilde{V}_1(z)$','interpreter','latex')
ylabel('$\tilde{V}_2(z)$','interpreter','latex')
zlabel('$\tilde{V}_3(z)$','interpreter','latex')
xlim([-2.5, 2.5]);
ylim([-2.5, 2.5]);
zlim([-2, 2.5]);
view([45 30]);
set_fonts(gcf, gca, fontsize);

%% RP using z reconstruction
subplot(3,4,4+7)
rp = crp(zz, zz, T, tau,radius,...
    'euclidean','nonormalize','silent');
rqa_z = crqa(zz, zz, T, tau,radius,...
    'euclidean','nonormalize','silent');
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set_fonts(gcf, gca, fontsize);


% MDRQA
subplot(3,4,4+8)
data = [xx, yy, zz];
%(DATA,DIM,EMB,DEL,NORM,RAD,ZSCORE)
[rp, mdrqa] = MDRQA(data, 3, 3, T, 'euc', 0.8*radius, 0);
rp = flipud(rp);
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set_fonts(gcf, gca, fontsize);

