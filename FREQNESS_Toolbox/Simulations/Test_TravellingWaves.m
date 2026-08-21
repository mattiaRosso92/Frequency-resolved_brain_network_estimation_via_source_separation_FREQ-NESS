ncomps =30;
% TIME, in seconds
t1 = 70;
t2 = 75;
srate = 250;

load("MNI152_8mm_coord_dyi.mat")

figure(9), clf
hold on
for compi = 1:ncomps

    % By default, expects 1 frequency in FREQ
    this_ts = FREQ{1}.ts(compi,srate*t1:srate*t2,1);

    filt_ts    = filterFGx(this_ts,srate,8.4,1.8,0);
    this_phase = angle(hilbert(filt_ts));

    this_pat = abs( FREQ{1}.pats(:,compi,1) );

    this_ampl = abs(hilbert(filt_ts));

    subplot(221)
    plot(this_ts)
    subplot(222)
    plot(filt_ts)
    hold on
    plot(this_ampl)
    hold off
    subplot(223)
    plot(this_phase)
    subplot(224)
    scatter3(MNI(:,1), MNI(:,2), MNI(:,3), 20, this_pat, 'filled');
    view(45, 20)
    colorbar
    colormap(parula)
    %clim([min(FREQ{1}.pats(:)), max(FREQ{1}.pats(:))])


    pause
end