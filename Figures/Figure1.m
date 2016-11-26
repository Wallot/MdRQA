%% Figure 1
% This script reproduces Figure 1 in:
%
% Wallot, S., Roepstorff, A., & Mønster, D. (2016). Multidimensional
% Recurrence Quantification Analysis (MdRQA) for the analysis of
% multidimensional time-series: A software implementation in MATLAB and its
% application to group-level data in joint action. Frontiers in Psychology,
% 7, 1835. http://dx.doi.org/10.3389/fpsyg.2016.01835
%
% Depends on:
% -Cross Recurrence Plot Toolbox
% -MdRQA script
%
% By:
% Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany
% Dan Mønster, Aarhus University, Aarhus, Denmark
%
% July 2016

fontsize = 15;

figure('units','normalized','position',[.1 .1 0.95 .35])
data=zscore(sin(0.1:0.1:20)+rand(1,200)/2);

subplot(2,3,1)
plot(data,'k')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\tilde{V}_1(t) = x(t)$', 'interpreter', 'latex')
set_fonts(gcf, gca, fontsize);

subplot(2,3,4)
plot(data(25:end),'k')
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\tilde{V}_2(t)$', 'interpreter', 'latex')
set_fonts(gcf, gca, fontsize);
axis([0 200 -2 2])

subplot(2,3,[2 5])
plot(data(1:end-16),data(17:end),'ok')
ax = gca;
xlabel('$\tilde{V}_1$', 'interpreter', 'latex')
ylabel('$\tilde{V}_2$', 'interpreter', 'latex')
ax.XTick = [-2 -1 0 1 2];
ax.YTick = [-2 -1 0 1 2];
axis square
set_fonts(gcf, gca, fontsize);

X=crp(data,2,15,.3);
subplot(2,3,[3 6])
spy(flipud(X), 'k')
ax = gca;
ax.XTick = [0 50 100 150];
ax.YTick = [0 50 100 150];
xlabel('')
set_fonts(gcf, gca, fontsize);
axis on
