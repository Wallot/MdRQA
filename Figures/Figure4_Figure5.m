% This script reproduces Figure 4 and Figure 5 in:
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
% -coupled_vdp.m
% -set_fonts.m
%
% Description:
% CRQA and MDRQA analysis of two coupled van der Pol oscillators
% This script applies Cross-Recurrence Quantification Analysis (CRQA) and
% Multi-Dimensional Recurrence Quantification Analysis (MDRQA) to a system
% of two coupled van der Pol oscillators.
%
% By:
% Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany
% Dan Mønster, Aarhus University, Aarhus, Denmark
%
% July 2016.
%

% Set the font size to use in plots.
fontsize = 12; 

crqa_zscore = true;
if crqa_zscore
    crqa_radius = 0.01;
else
    crqa_radius = 0.01;
end
mdrqa_radius = 0.01;
crqa_embed.m = 2;
crqa_embed.t = 1;
mdrqa_embed.m = 1;
mdrqa_embed.t = 1; % Only used for mdrqa_embed.m > 1

x0 = 1;
y0 = 2;
xdot0 = 0;
ydot0 = 0;
mu = 100;
e1 = 0.01;
e2 = 0.05;
tmax = 5.e3;


[t x y] = coupled_vdp([x0, y0, xdot0, ydot0], mu, e1, e2, tmax);

%% Plot the time series
figure('units','normalized','position',[.1 .1 .75 .6])
subplot(2,2,1)
plot(x,'MarkerFaceColor','b','linewidth',1, 'markers',4);
hold on
plot(y,'MarkerFaceColor','r','linewidth',1, 'markers',4);
xlabel('$t$','interpreter','latex');
ylabel('$x(t), y(t)$','interpreter','latex');
set_fonts(gcf, gca, fontsize);


%% Resample the time series
Nsamples = 1e3;
sampling_interval = tmax/Nsamples;
[tx ty] = synchronize(x,y,'Uniform','Interval',sampling_interval);

%% Low coupling strength
% CRQA
subplot(2,4,5)
if crqa_zscore
    xx = zscore(tx.Data);
    yy = zscore(ty.Data);
else
    xx = tx.Data;
    yy = ty.Data;
end
rp = crp(xx, yy,crqa_embed.m,crqa_embed.t,crqa_radius,...
    'euclidean','nonormalize','silent');
crqa_low_coupling = crqa(xx, yy, crqa_embed.m, crqa_embed.t, crqa_radius,...
    'euclidean','nonormalize','silent');
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set(gca,'FontSize',fontsize,'fontWeight','normal')
set(findall(gcf,'type','text'),'FontSize',fontsize,'fontWeight','normal')

% MDRQA
data = [xx,yy];
[rp, mdrqa_low_coupling] = MDRQA(data,2,...
    mdrqa_embed.m, mdrqa_embed.t,'euc', mdrqa_radius, 0, 0);
subplot(2,4,6)
rp = flipud(rp);
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set_fonts(gcf, gca, fontsize);


%% High coupling strength
x0 = 2;
y0 = 1;
xdot0 = 0;
ydot0 = 0;
mu = 100;
e1 = 0.02;
e2 = 0.1;
tmax = 5.e3;
% Generate new data
[t x y] = coupled_vdp([x0, y0, xdot0, ydot0], mu, e1, e2, tmax);

% Plot the time series
subplot(2,2,2)

plot(x,'MarkerFaceColor','b','linewidth',1, 'markers',4);
hold on
plot(y,'MarkerFaceColor','r','linewidth',1, 'markers',4);
%pbaspect([1.6 1 1]);
xlabel('$t$','interpreter','latex');
ylabel('$x(t), y(t)$','interpreter','latex');
set_fonts(gcf, gca, fontsize);

%% temp section
% Resample the time series
[tx ty] = synchronize(x,y,'Uniform','Interval',sampling_interval);
if crqa_zscore
    xx = zscore(tx.Data);
    yy = zscore(ty.Data);
else
    xx = tx.Data;
    yy = ty.Data;
end

% CRQA
subplot(2,4,7)
rp = crp(xx, yy,crqa_embed.m,crqa_embed.t,crqa_radius,...
    'euclidean','nonormalize','silent');
crqa_high_coupling = crqa(xx, yy, crqa_embed.m, crqa_embed.t, crqa_radius,...
    'euclidean','nonormalize','silent');
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set_fonts(gcf, gca, fontsize);

% MDRQA
data = [xx,yy];
[rp, mdrqa_high_coupling] = MDRQA(data,2,...
    mdrqa_embed.m, mdrqa_embed.t, 'euc', mdrqa_radius, 0, 0);
subplot(2,4,8)
rp = flipud(rp);
spy(rp, 'k');
axis xy;
ax = gca;
ax.XTick = [0 500 1000];
ax.YTick = [0 500 1000];
xlabel('')
set_fonts(gcf, gca, fontsize);


%% plot RP measures (e.g., ADL) as a function of coupling strength.
epsilon_values = [0.01:0.02:0.2]';
if exist('rqa_results.mat', 'file') ~= 2
    N_sample = 1;
    crqa_results = zeros(numel(epsilon_values), numel(crqa_low_coupling));
    mdrqa_results = zeros(numel(epsilon_values), numel(mdrqa_low_coupling));
    for s = 1:N_sample
        for k = 1:numel(epsilon_values)
            e2 = epsilon_values(k);
            e1 = e2/5.0;
            % Calculate timeseries
            [t x y] = coupled_vdp([x0, y0, xdot0, ydot0], mu, e1, e2, tmax);
            % Resample the time series
            [tx ty] = synchronize(x,y,'Uniform','Interval',sampling_interval);
            if crqa_zscore
                xx = zscore(tx.Data);
                yy = zscore(ty.Data);
            else
                xx = tx.Data;
                yy = ty.Data;
            end
            % CRQA
            crqa_result = crqa(xx, yy, crqa_embed.m, crqa_embed.t,...
                crqa_radius,'euclidean','nonormalize','silent');
            crqa_results(k,:) = crqa_results(k,:) + crqa_result;
            % MDRQA
            data = [xx, yy];
            [rp, mdrqa_result] = MDRQA(data,2, mdrqa_embed.m,...
                mdrqa_embed.t, 'euc', mdrqa_radius, 0, 0);
            mdrqa_results(k,:) = mdrqa_results(k,:) + mdrqa_result;
        end
    end
    crqa_results = crqa_results / N_sample;
    mdrqa_results = mdrqa_results / N_sample;
    save('rqa_results.mat', 'crqa_results', 'mdrqa_results')
else
    load('rqa_results.mat');
end

%% Plot results
figure()
rqa_measures = containers.Map({'RR (%)', 'DET (%)', 'ADL', 'LDL', 'ENT', 'LAM',...
    'TT','LVL', 'T1', 'T2'},...
    {1,2,3,4,5,6,7,8,9,10});
measures = {'RR (%)', 'DET (%)', 'ADL', 'LDL'};
eps1 = epsilon_values/5.0;
for k = 1:numel(measures)
    subplot(2,2,k)
    measure = measures{k};
    disp(sprintf('Results for %s:', measure))
    idx = rqa_measures(measure);
    % Factor to convert RR to percent
    if (k == 1 | k == 2)
        factor = 100.0;
    else
        factor = 1.0;
    end
    plot(eps1(:),crqa_results(:,idx)*factor,'--k', 'linewidth',1)
    hold on
    % MdRQA measures offset by 1 relative to CRQA, so idx + 1 is used.
    plot(eps1(:),mdrqa_results(:,idx+1),'-k', 'linewidth',1) 
    % Make sure 0 is included on y-axis
    yrange = ylim();
    if k == 2
        ylim([80 yrange(2)]);
    else
        ylim([0 yrange(2)]);
    end
   
    xlabel('$\epsilon_1$','interpreter','latex');
    ylabel(measure);
    disp('Corr(CRQA, MdRQA):')
    rho = corr(crqa_results(:,idx),mdrqa_results(:,idx+1),...
        'type', 'Pearson');
    disp(rho)
    disp('Corr(CRQA, eps1):')
    rho = corr(crqa_results(:,idx),eps1(:),...
        'type', 'Pearson');
    disp(rho)
    disp('Corr(MdRQA, eps1):')
    rho = corr(mdrqa_results(:,idx+1),eps1(:),...
        'type', 'Pearson');
    disp(rho)
    label_string = sprintf('$\\rho = %.2g$', rho);
    xrange = xlim();
    yrange = ylim();
    set_fonts(gcf, gca, fontsize);
end


