function analyze_MSD_directed_motion(vertex_tracks)

SPACE_UNITS = 'µm';
TIME_UNITS = 's';

tracks = vertex_tracks

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

ma = ma.addAll(tracks);
% 
% figure
% ma.plotTracks;
% ma.labelPlotTracks;
% 
ma = ma.computeMSD;
ma.msd
% 
% figure
% ma.plotMSD;
% 
% figure
% cla
% ma.plotMeanMSD(gca, true)

A = ma.getMeanMSD;
t = A(:, 1); % delay vector
msd = A(:,2); % msd
std_msd = A(:,3); % we will use inverse of the std as weights for the fit
std_msd(1) = std_msd(2); % avoid infinity weight

ft = fittype('a*x + c*x^2');
[fo, gof] = fit(t, msd, ft, 'Weights', 1./std_msd, 'StartPoint', [0 0]);

hold on
plot(fo)
legend off
ma.labelPlotMSD

Dfit = fo.a / 4;
Vfit = sqrt(fo.c);

ci = confint(fo);
Dci = ci(:,1) / 4;
Vci = sqrt(ci(:,2));

fprintf('Parabolic fit of the average MSD curve with 95%% confidence interval:\n')

fprintf('D = %.3g [ %.3g - %.3g ] %s\n', ...
    Dfit, Dci(1), Dci(2), [SPACE_UNITS '²/' TIME_UNITS]);

fprintf('V = %.3g [ %.3g - %.3g ] %s\n', ...
    Vfit, Vci(1), Vci(2), [SPACE_UNITS '/' TIME_UNITS]);


%%%%%

ma = ma.computeDrift('velocity');
% figure
% ma.plotDrift
% ma.labelPlotTracks

ma = ma.computeMSD;
figure
ma.plotMeanMSD(gca, true)
[fo gof] = ma.fitMeanMSD(0.75);
plot(fo)
ma.labelPlotMSD;
legend off

end
