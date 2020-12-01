% Tim C Whalen, last edited July 2020
% Makes Whalen 2020 JNeurophys Figure 7, GPe and STN fraction of delta oscillating units

allregions = {'GPe' 'STN'};
allconds = {'control' 'acute'};
conds_full = {'Control' 'Bilateral'};
nconds = length(allconds);
nregions = length(allregions);
fracs = cell(nconds,nregions);

for c = 1:length(allconds)
    for r = 1:length(allregions)
        load(['Whalen2020_data_' allregions{r} '_' allconds{c}],'data');
        data = renewalPSD_phaseShift_batch(data);
        fracs{c,r} = data.osc.frac_osc;
    end
end

% colorblind-friendly palette from Wong 2011 Nature Methods
cols = [76 76 76;
    0 158 115; % teal
    86 180 233; % sky blue
    230 159 0; % orange
    0 114 178; % blue
    123 94 0]./255; % vermillion
cols_dark = cols./1.3;
cols_dark(1,:) = [0 0 0];
shapes = 'od^s';

figure
subplot(nregions,1,1)
hold on
for r = 1:nregions
    subplot(nregions,1,r)
    hold on
    for c = 1:nconds
        ps = fracs{c,r};
        np = length(ps);
        me = mean(ps);
        plot([-.1 .1]+c,[me me],'Color',1-(1-cols(c,:))/1.5,'LineWidth',2)
        scatter(c+((0:.5:(np-1)/2)-np/4)/10, ps, 15,cols(c,:),shapes(c),'MarkerEdgeColor',cols_dark(c,:),'MarkerFaceColor',cols(c,:))
    end
    ylim([0 1])
    xlim([.5 4.5])
    if r == nregions
        set(gca, 'xtick',1:nconds,'XTickLabels',{'Naive' 'Bilateral' 'Unilateral'},'ytick', [0 1]) %,'xcolor','w')
    else
        set(gca, 'xtick',[],'ytick', [0 1])
    end
    title([allregions{r} ' Power + Phase Shift (0.5-4 Hz)'])
    makeNice(gca)
end

figure
subplot(2,1,1)
hold on
plotFracs(fracs(:,1),cols,shapes,'GPe Power + Phase Shift (0.5-4 Hz)',conds_full)
subplot(2,1,2)
hold on
plotFracs(fracs(:,2),cols,shapes,'STN Power + Phase Shift (0.5-4 Hz)',conds_full)

[~, p_gpe] = ttest2(fracs{1,1},fracs{2,1})
[~, p_stn] = ttest2(fracs{1,2},fracs{2,2})