% Tim C Whalen, last edited July 2020
% Makes Whalen 2020 JNeurophys Figure 3, fraction of delta oscillating units
allregions = {'snr'};
allconds = {'control' 'acute' 'uni_ipsi' 'reserpine'};
conds_full = {'Control' 'Bilateral' 'Unilateral' 'Reserpine'};
nconds = length(allconds);
nregions = length(allregions);
fracs = cell(nconds,nregions);
fracs_old = cell(nconds,nregions);
fracs_beta = cell(nconds,nregions);
fracs_old_beta = cell(nconds,nregions);

for c = 1:length(allconds)
    for r = 1:length(allregions)
        load(['Whalen2020_data_SNr_' allconds{c}],'data');
        data.osc = struct();
        data_beta = data; % compute beta separately
        data_beta.osc.srch_lo = 7; % note delta is function's default
        data_beta.osc.srch_hi = 35;
        
        data = renewalPSD_phaseShift_batch( data);
        fracs{c,r} = data.osc.frac_osc;
        fracs_old{c,r} = data.osc.frac_osc_old;
        data_beta = renewalPSD_phaseShift_batch( data_beta);
        fracs_beta{c,r} = data_beta.osc.frac_osc;
        fracs_old_beta{c,r} = data_beta.osc.frac_osc_old;
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
shapes = 'od^sv';

figure
subplot(2,2,1)
hold on
plotFracs(fracs,cols,shapes,'Power + Phase Shift (0.5-4 Hz)',conds_full)
subplot(2,2,3)
hold on
plotFracs(fracs_old,cols,shapes,'Power Only (0.5-4 Hz)',conds_full)
subplot(2,2,2)
hold on
plotFracs(fracs_beta,cols,shapes,'Power + Phase Shift (7 -35 Hz)',conds_full)
subplot(2,2,4)
hold on
plotFracs(fracs_old_beta,cols,shapes,'Power Only (7 -35 Hz)',conds_full)

groups = cell(size(fracs));
for i = 1:size(fracs,1)
    groups{i} = i+zeros(size(fracs{i}));
end
groups = cell2mat(groups);

% tests - param
matfracs = cell2mat(fracs);
[pan, ~, stats] = anova1(matfracs,groups,0)
if pan < 0.05
    ps = dunnett(stats,2:4,1);
end

matfracs_old = cell2mat(fracs_old);
[pan_old, ~, stats_old] = anova1(matfracs_old,groups,0)
if pan_old < 0.05
    ps_old = dunnett(stats_old,2:4,1);
end

matfracs_beta = cell2mat(fracs_beta);
[pan_beta, ~, stats_beta] = anova1(matfracs_beta,groups,0)
if pan_beta < 0.05
    ps_beta = dunnett(stats_beta,2:4,1);
end

matfracs_old_beta = cell2mat(fracs_old_beta);
[pan_old_beta, ~, stats_old_beta] = anova1(matfracs_old_beta,groups,0)
if pan_old_beta < 0.05
    ps_old_beta = dunnett(stats_old_beta,2:4,1);
end