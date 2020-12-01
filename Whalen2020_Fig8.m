% Tim C Whalen, last edited Nov 2020
% Makes Whalen et al. 2020 JNeurophys Figure 8, M1 ECoG analysis

%% Panels B-C - M1 ECoG bandpower
examples = [19; 73];
wind=2^14; nover=2^13;
maxtime = 360; % baseline time, cuts off any experimental manipulations
ncond = 2;

% load only ecog (regardless of deep location, as no spikes are included)
structs = cell(ncond,1);
load('Whalen2020_data_ecog_control','data');
structs{1} = data;
load('Whalen2020_data_ecog_acute','data');
structs{2} = data;
clear data

fractheta = cell(ncond,1);
fracdelta = cell(ncond,1);

figure
for i = 1:ncond
    data = structs{i};
    FS = data.ECOG_FS(1);

    fracdelt = zeros(data.nfiles,1);
    fracthet = zeros(data.nfiles,1);
    for f = 1:data.nfiles
        ecog = data.ecogs_nofilt{f};
        powdelt = bandpower(ecog,data.ECOG_FS(f),[.5 4]);
        powthet = bandpower(ecog,data.ECOG_FS(f),[4 8]);
        powfull = bandpower(ecog,data.ECOG_FS(f),[.5 100]);
        fracdelt(f) = powdelt/powfull;
        fracthet(f) = powthet/powfull;
    end
    fractheta{i} = zeros(max(data.animalcodes),1);
    fracdelta{i} = zeros(max(data.animalcodes),1);
    for a = 1:max(data.animalcodes)
        files = find(data.animalcodes==a);
        fractheta{i}(a) = mean(fracthet(files));
        fracdelta{i}(a) = mean(fracdelt(files));
    end
    
    map = 2*(i-1)+1; % quick hack for subplotting
    subplot(1,5,map:map+1)
    f = examples(i);
    ecog = data.ecogs_nofilt{f};
    if length(ecog) > maxtime*FS+1
        ecog = ecog(1:maxtime*FS+1);
        T = maxtime;
    else
        T = data.T(f);
    end
    [pow, freqs] = pwelch(ecog,wind,nover,[],FS);
    plot(freqs,pow/T/bandpower(ecog,data.ECOG_FS(f),[.5 100]),'k','LineWidth',2)
    xlim([0 20])
    ylim([0 2*10^-3])
    makeNice(gca)
end

barmat = {fracdelta{1} fracdelta{2}; fractheta{1} fractheta{2}};
subplot(1,5,5)
barwitherr(cellfun(@(x) std(x)/sqrt(length(x)), barmat),cellfun(@mean,barmat))
ylim([0 1])
makeNice(gca)

%% Panels E-F - example AP/IP regression results
load('Whalen2020_data_SNr_acute_ecog.mat','data')
example_files = [21,27];
data_ex = makeSubstruct(data, example_files); % only process needed files from DD data
data_ex.ecog_reg = struct('nlag_past',100,'nlag_fut',100,'ar_lag',10);
[ data_ex ] = ecogSpikeRegress_batch(data_ex);

nlag_past = data_ex.ecog_reg.nlag_past;
nlag_fut = data_ex.ecog_reg.nlag_fut;
nlag_step = data_ex.ecog_reg.step;
step = data_ex.ecog_reg.step;
bs_all = data_ex.ecog_reg.bs;
ssrs_all = data_ex.ecog_reg.ssrs;

f = 1;  u = 2; % example for AP
figure
subplot(2,8,1:3)
hold on
plot(1000*(-nlag_past*step:step:nlag_fut*step),bs_all{f}(u,:),'k','LineWidth',2)
ylabel('Regression Coeffs.')
xlabel('<-- Past Spikes           Lag (msec)         Future Spikes -->')
makeNice(gca)

subplot(2,8,9:11)
plot(1000*(-nlag_past*step:step:nlag_fut*step),ssrs_all{f}(u,:),'k','LineWidth',2)
ylabel('\Sigma (residuals^2)')
lag = data_ex.ecog_reg.sig_lags{f}(u);
if ~isnan(lag)
    hold on
    plot(1000*lag, min(ssrs_all{f}(u,:)), 'b.','MarkerSize',15)
    text(1000*(lag+.02),min(ssrs_all{f}(u,:)),int2str(round(lag*1000)),'Color','r','FontSize',10)
end
makeNice(gca)

f = 2;  u = 2; % example for IP
subplot(2,8,4:6)
hold on
plot(1000*100*(-nlag_past*step:step:nlag_fut*step),bs_all{f}(u,:),'k','LineWidth',2)
makeNice(gca)

subplot(2,8,12:14)
plot(1000*(-nlag_past*step:step:nlag_fut*step),ssrs_all{f}(u,:),'k','LineWidth',2)
lag = data_ex.ecog_reg.sig_lags{f}(u);
if ~isnan(lag)
    hold on
    plot(1000*lag, min(ssrs_all{f}(u,:)), 'r.','MarkerSize',15)
    text(1000*(lag+.02),min(ssrs_all{f}(u,:)),int2str(round(lag*1000)),'Color','r','FontSize',10)
end
makeNice(gca)

%% Figure 8G-H - histograms of coefficients and lags

data.osc=struct('min_rate', 0);
[data] = renewalPSD_phaseShift_batch(data);
[data] = ecogSpikeRegress_batch(data);
pttypebu = zeros(1,0); % bu = by unit. 1 if peak first, -1 if trough
troughbu = zeros(1,0);
peakbu = zeros(1,0);
sig_lagsbu = zeros(1,0);
sig_freqsbu = zeros(1,0);
coef_at_sig = zeros(1,0);
anbu = zeros(1,0); % animal
count = 0;
for f = 1:data.nfiles
    nu = size(data.ts{f},1);
    for u = 1:nu
        if data.osc.has_osc{f}(u)
            count = count+1;
            pttypebu(count,1) = data.ecog_reg.pttype{f}(u);
            sig_lagsbu(count,1) = data.ecog_reg.sig_lags{f}(u);
            anbu(count,1) = data.animalcodes(f);
            sig_freqsbu(count,1) = data.osc.freqs(data.osc.sig_inds{f}{u}(1)); % assumes only one sig freq found
            if isnan(data.ecog_reg.sig_lags{f}(u))
                coef_at_sig(count,1) = NaN;
            else
                coef_at_sig(count,1) = data.ecog_reg.bs{f}(u,data.ecog_reg.nlag_past+1+round(data.ecog_reg.sig_lags{f}(u)/data.ecog_reg.bin))/data.ecog_reg.mean_amp(f); % converts lag in msec to index of bs, noting that lag may be negative
            end
        end
    end
end
sig_phasebu = 2*pi*sig_lagsbu.*sig_freqsbu;

edges = [-3.5:.25:3.5]*10^-4;
bpeak = histcounts(coef_at_sig(pttypebu==1),edges);
btrough = histcounts(coef_at_sig(pttypebu==-1),edges);

subplot(2,8,7:8)
hold on
h = bar(edges(1:end-1),bpeak,'histc');
h.FaceColor = [.8 0 0];
h.LineStyle = 'none';
h.FaceAlpha = 0.5;
h = bar(edges(1:end-1),btrough,'histc');
h.FaceColor = [0 0 .8]*0.8;
h.LineStyle = 'none';
h.FaceAlpha = 0.5;
plot([0 0],[0 10],'k--','LineWidth',2)
legend('Peak','Trough')
xlabel('Regression Coeff.')
ylabel('Unit Count')
xlim([min(edges) max(edges)])
makeNice(gca)

edges = -pi:pi/16:0;
peakbar = zeros(max(data.animalcodes),length(edges)-1);
troughbar = zeros(max(data.animalcodes),length(edges)-1);

subplot(2,8,15:16)
peakbar = histcounts(sig_phasebu(pttypebu==1  & sig_lagsbu<0),edges);
troughbar = histcounts(sig_phasebu(pttypebu==-1 & sig_lagsbu<0),edges);
hold on
h = bar(edges(1:end-1),peakbar,'histc');
h.FaceColor = [.8 0 0];
h.LineStyle = 'none';
h.FaceAlpha = 0.5;
set(gca,'XTick',[-pi:pi/2:0],'XTickLabels',{'\pi' '\pi/2' '0'})
hold on
h = bar(edges(1:end-1),troughbar,'histc');
h.FaceColor = [0 0 0.8];
h.LineStyle = 'none';
h.FaceAlpha = 0.5;
set(gca,'XTick',[-pi:pi/2:0],'XTickLabels',{'\pi' '\pi/2' '0'})
xlim([-pi pi/2])
xlabel('Phase Offset from M1 (rad)')
ylabel('Unit Count')
makeNice(gca)

%% Panel I - AP/IP rate and CV
ts = getit(data.ts)';
ts_osc = ts(cell2mat(data.osc.has_osc)==1);
rates = cell2mat(data.rates);
rates_osc = rates(cell2mat(data.osc.has_osc)==1);

minr = 5;
ts_ap = ts_osc(pttypebu==1 & rates_osc>minr);
ts_ip = ts_osc(pttypebu==-1 & rates_osc>minr);

ratefun = @(x) length(x)/(max(x)-min(x));
rates_ap = cellfun(ratefun,ts_ap);
rates_ip = cellfun(ratefun,ts_ip);

cvfun = @(x) std(diff(x))/mean(diff(x));
cv_ap = cellfun(cvfun,ts_ap);
cv_ip = cellfun(cvfun,ts_ip);

ISI_ap = cellfun(@diff, ts_ap, 'UniformOutput',0);
ISI_ip = cellfun(@diff, ts_ip, 'UniformOutput',0);

figure
subplot(1,2,1)
boxplot([rates_ap;rates_ip],[zeros(size(cv_ap));ones(size(cv_ip))],'symbol','')
title('Rates')
ylabel('Firing Rate (Hz)')
ylim([0 40])
set(gca,'xtick',[1 2],'xticklabels',{'AP' 'IP'})
makeNice(gca)
subplot(1,2,2)
boxplot([cv_ap;cv_ip],[zeros(size(cv_ap));ones(size(cv_ip))],'symbol','')
title('CV ISI')
ylabel('CV ISI')
ylim([0 4])
set(gca,'xtick',[1 2],'xticklabels',{'AP' 'IP'})
makeNice(gca)