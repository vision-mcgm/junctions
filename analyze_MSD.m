function analyze_MSD(vertex_tracks)

SPACE_UNITS = '�m';
TIME_UNITS = 's';

tracks = vertex_tracks;

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

ma = ma.addAll(tracks);

figure
ma.plotTracks;
ma.labelPlotTracks;

ma = ma.computeMSD;
ma.msd

figure
ma.plotMSD;

figure
cla
ma.plotMeanMSD(gca, true)

% mmsd = ma.getMeanMSD;
% t = mmsd(:,1);
% x = mmsd(:,2);
% dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
% errorbar(t, x, dx, 'k')

[fo, gof] = ma.fitMeanMSD(0.75);
plot(fo)
ma.labelPlotMSD;
legend off

%fitting each individual curve
ma = ma.fitMSD(0.75);



end
