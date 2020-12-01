% Tim C Whalen, last edited July 2020
% Makes Whalen 2020 JNeurophys Figure 9, oscillations with M1 lesion

allregions = {'snr'};
allconds = {'control' 'acute' 'asp'};
allexamples = {[] [] [4 1]}; % [file, unit] examples to plot
conds_full = {'Control' 'DD' 'DD + M1 Aspiration'};
nconds = length(allconds);
nregions = length(allregions);
fracs = cell(nconds,nregions);

for c = 1:length(allconds)
    for r = 1:length(allregions)
        load(['Whalen2020_data_SNr_' allconds{c}],'data');
        data.osc = struct('examples',allexamples{c});
        [ data ] = renewalPSD_phaseShift_batch( data);
        fracs{c,r} = data.osc.frac_osc;
    end
end

cols = [76 76 76;
    0 158 115; % teal
    123 94 0]./255; % vermillion
cols_dark = cols./1.3;
cols_dark(1,:) = [0 0 0];
shapes = 'od^';

figure
hold on
plotFracs(fracs,cols,shapes,'Power + Phase Shift (0.5-4 Hz)',conds_full)

groups = cell(size(fracs));
for i = 1:size(fracs,1)
    groups{i} = i+zeros(size(fracs{i}));
end
groups = cell2mat(groups);

% tests
matfracs = cell2mat(fracs);
[pan, ~, stats] = anova1(matfracs,groups,0)
if pan < 0.05
    ps = dunnett(stats,1:2,3);
end

[h, p] = ttest2(fracs{1},fracs{2})