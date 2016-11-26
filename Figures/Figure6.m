% This script reproduces Figure 6 in:
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
% This script applies Joint Recurrence Quantification Analysis (JRQA) and
% Multi-Dimensional Recurrence Quantification Analysis (MDRQA) to the
% Lorenz model.
%
% By:
% Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany
% Dan Mønster, Aarhus University, Aarhus, Denmark
%
% September 2016.
%

% Parameters:
fontsize = 14; % Font size to use in plots.
x0 = 10;
y0 = 10;
z0 = 10;
sigma = 10;
rho = 28;
beta = 8/3.0;
tmax = 20;

[t, x, y, z] = lorenz_model([x0, y0 z0], sigma, rho, beta, tmax);
xx = zscore(x.Data);
yy = zscore(y.Data);
zz = zscore(z.Data);



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
E = 3;
tau = 4;
% RQA radius
rqa_radius = 0.1;
mdrqa_radius = 0.8 * rqa_radius;

%% RP using x, y, and z reconstructions

rp_x = crp(xx, xx, E, tau,rqa_radius,'euclidean','nonormalize','silent');
rp_y = crp(yy, yy, E, tau,rqa_radius,'euclidean','nonormalize','silent');
rp_z = crp(zz, zz, E, tau,rqa_radius,'euclidean','nonormalize','silent');

% MDRQA
data = [xx, yy, zz];
[rp, mdrqa] = MDRQA(data,3, E, 0, 'euc', mdrqa_radius, 0, 0);
rp = flipud(rp);

%%  Plot JRP and MdRP
%Plot the joint recurrence plot obtained from all three variables x, y, and
% z; and the multidimensional recurrence plot.
figure()
subplot(1,2,1)
jrp=rp_x.*rp_y.*rp_z;
spy(jrp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
tit = title('a)');
tit.Position = [0 1300 0];
set_fonts(gcf, gca, fontsize);

subplot(1,2,2)
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
tit = title('b)');
tit.Position = [0 1300 0];
set_fonts(gcf, gca, fontsize);
