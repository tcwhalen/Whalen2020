% Tim C Whalen, last edited Nov 2020
% Makes Whalen 2020 JNeurophys Figure 5 summary figures, dynamic
% relationship between oscillation power and motor activity

load('Whalen2020_data_fig5')

% quick fix to scale new mv to amplitude of old mv
data.mv_proc = cellfun(@(x) x./.0192,data.mv_proc,'UniformOutput',0);

% get baseline oscillating units
data.osc = struct();
data.osc.wind=2^13;
data = renewalPSD_phaseShift_batch( data);
osc_only = 1; % only units with detected oscillation across whole time

try_maxlag = 50; % for xcorr between spectrogram and movement, "try" as it will be rounded to nearest multiple of step, sec
psd_thresh = .1; % threshold for finding index of osc peak, only used if osc_only==0
preproc_thresh = 5; % threshold for file containing enough movement to be used
ismv_thresh = .02; % mvmt threshold on mean over several bins for defining bouts (rest must be exactly 0)
NFFT = data.osc.wind;
step = NFFT/4;
srch_lo = .5;
srch_hi = 4;

FS = 1000; % to downsample spike trains to
freqs = 0:FS/NFFT:srch_hi+1;
srch_inds_full = find(freqs>= srch_lo & freqs <= srch_hi);

max_pow_rest = []; % max 0.5-4 Hz power during rest, by unit
max_pow_mvmt = [];
nseg_rest = []; % # segments used to calculate rest PSD
nseg_mvmt = [];
max_freq_rest = [];
max_freq_mvmt = [];
psds_all = cell(data.nfiles,1);
oscmov_all = cell(data.nfiles,1);
mmeans_all = cell(1,0);

for f = 1:data.nfiles
    if ~isempty(data.movet_raw{f})
        movet = data.movet_raw{f}-data.movet_rs{f}(1);
        movey = data.movey_raw{f};
    else
        movet = data.mv_proctime{f};
        movey = data.mv_proc{f};
    end
    mvend = movet(end);
    
    if sum(movey>1) < preproc_thresh
        psds_all{f} = NaN;
        oscmov_all{f} = NaN;
        continue
    end
    if osc_only
        ts_all = data.ts{f};
        ts_rate = ts_all(data.rates{f}>data.osc.min_rate);
        ts = ts_rate(data.osc.has_osc{f}==1);
    else
        ts = data.ts{f};
    end
    sig_inds_all = data.osc.sigp_inds{f}(data.osc.has_osc{f}==1);
    
    nu = length(ts);
    psds_all{f} = cell(nu,1);
    oscmov_all{f} = cell(nu,1);
    for u = 1:length(ts)
        spk_t = ts{u};
        spk_t = spk_t(spk_t<mvend);
        
        spk.t = round( FS.*spk_t );
        len_delt = spk.t(end)-spk.t(1)+1;
        spk.delt = zeros(1,len_delt);
        spk.delt( spk.t - spk.t(1)+1 ) = 1;
        
        Ncalc = floor(len_delt/step - NFFT/step +1); % # fft's to do
        psds = zeros(NFFT/2,Ncalc);
        full_freqs = 0:FS/NFFT:FS/2;
        min_norm = find(full_freqs>=200);
        windfun = hamming(NFFT)';
        mmeans = zeros(1,Ncalc); % movement means
        
        for s = 1:Ncalc
            mseg = abs(movey(movet>=(1+step*(s-1))/FS & movet<(step*(s-1)+NFFT)/FS));
            mseg(mseg==1)=0;
            mmeans(s) = mean(mseg);
            delt = spk.delt(1+step*(s-1):step*(s-1)+NFFT);
            % Compare to renewal process defined by ISIs
            isi = diff(find(delt));
            edges = (0:NFFT)-.001;
            isin = histcounts(isi,edges);
            isidist = isin/sum(isin);
            phat = fft(isidist);
            chat = real((1+phat)./(1-phat));
            ff = abs(fft(windfun.*(delt-mean(delt)))).^2/FS/2;
            psd_corr = ff./chat;
            psd_corr = psd_corr./mean(psd_corr(min_norm:NFFT/2)); % normalize trick

            if sum(isnan(psd_corr))>0
                psd_corr = ones(size(psd_corr));
            end
            psds(:,s) = psd_corr(1:NFFT/2);
        end
        psds_all{f}{u} = psds;
        
        mz = zscore(mmeans);
        mmeans_all{end+1} = mmeans;
        
        maxi = ceil(try_maxlag/(step/FS));
        maxlag = maxi*(step/FS);
        xcorrt = -maxlag:step/FS:maxlag;
        
        if osc_only
            srch_inds = [];
            max_change = 4;
            sig_inds = sig_inds_all{u};
            for i = 1:length(sig_inds)
                srch_inds = [srch_inds sig_inds(i)-max_change:sig_inds(i)+max_change];
            end
            srch_inds = unique(srch_inds);
        else
            srch_inds = srch_inds_full;
        end
        
        % find avg psd and peak for rest
        psd_avg_rest = mean(psds(:,mmeans==0),2);
        psd_avg_rest = conv(psd_avg_rest,[1 1]./2, 'same');
        try
            if osc_only
                [~,peak_rest_ii] = max(psd_avg_rest(srch_inds));
                peak_rest = srch_inds(peak_rest_ii);
            else
                peak_rest = find_sig_peaks(psd_avg_rest(1:srch_inds(end)+3), psd_thresh, 3, srch_inds);
            end
        catch % hits if no rest periods, which is fine
            continue;
        end
        [max_pow_rest(end+1),max_indind_rest] = max(psd_avg_rest(peak_rest)); % index of sig_inds which are themselves indices
        max_ind_rest = peak_rest(max_indind_rest); % may be multiple if osc_only==0, so find best
        max_freq_rest(end+1) = freqs(max_ind_rest);
        
        % same for mvmt
        psd_avg_mvmt = mean(psds(:,mmeans>ismv_thresh),2);
        psd_avg_mvmt = conv(psd_avg_mvmt,[1 1]./2, 'same');
        nseg_mvmt(end+1) = sum(mmeans>ismv_thresh);
        nseg_rest(end+1) = sum(mmeans==0);
        do_mvmt_psd = 0;
        if sum(isnan(psd_avg_mvmt))==0 % if there is a movement period > thresh
            do_mvmt_psd =1;
            if osc_only
                [~,peak_mvmt_ii] = max(psd_avg_mvmt(srch_inds));
                peak_mvmt = srch_inds(peak_mvmt_ii);
            else
                peak_mvmt = find_sig_peaks(psd_avg_mvmt(1:srch_inds(end)+3), psd_thresh, 3, srch_inds);
            end
            [max_pow_mvmt(end+1),max_indind_mvmt] = max(psd_avg_mvmt(peak_mvmt)); % index of sig_inds which are themselves indices
            max_ind_mvmt = peak_mvmt(max_indind_mvmt);
            max_freq_mvmt(end+1) = freqs(max_ind_mvmt);
        else
            max_pow_mvmt(end+1)=NaN;
            max_freq_mvmt(end+1)=NaN;
        end
       
        cents = (NFFT/2:step:NFFT/2+step*(Ncalc-1))/FS;
        
        lowinds_rpeak = [max_ind_rest-1 max_ind_rest+1];
        lowfreq = mean(psds(lowinds_rpeak(1):lowinds_rpeak(2),:),1);
        mp_xcorr_rpeak = xcorr(mz, zscore(lowfreq),maxlag/(step/FS));
        if do_mvmt_psd
            lowinds_mvpeak = [max_ind_mvmt-1 max_ind_mvmt+1];
            lowfreq = mean(psds(lowinds_mvpeak(1):lowinds_mvpeak(2),:),1);
            mp_xcorr_mpeak = xcorr(mz, zscore(lowfreq),maxlag/(step/FS));
        end
    end
end

bad_bin = zeros(size(nseg_rest));
bad_bin(bad_units)=1;
seg_min = 8;
good_units = (nseg_rest>=seg_min) & (nseg_mvmt>=seg_min) & ~bad_bin;

good_pow_rest = max_pow_rest(good_units==1);
good_pow_mvmt = max_pow_mvmt(good_units==1);
good_pow_diff = good_pow_mvmt-good_pow_rest;
edges = -6:.5:6;

figure
hold on
hi = histogram(good_pow_diff,edges);
hi.FaceColor = [0 0 .5];
hi.EdgeColor = [0 0 .5];
maxval = 23; % somewhere around max of histogram
dotscale = 5;

plot([0 0],[0 maxval+dotscale+2],'--','Color',[.8 .8 .8],'LineWidth',2)
scatter(good_pow_diff,maxval+dotscale*rand(size(good_pow_diff)),15,'.','b')
m = mean(good_pow_diff);
plot([m m],[maxval-1 maxval+dotscale+1],'k','LineWidth',2)
ylim([0 maxval+dotscale+1])
ylabel('\Delta Power (Movement-Rest)')
xlabel('Count')

figure
hold on
edges_pos = 0:1:17.5;
hir = histogram(good_pow_rest,edges_pos);
hir.FaceColor = [.1 .1 .1];
hir.EdgeColor = [.1 .1 .1];
hir.FaceAlpha = 0.5;
hir.EdgeAlpha = 0.5;
him = histogram(good_pow_mvmt,edges_pos);
him.FaceColor = [.7 0 0];
him.EdgeColor = [.7 0 0];
him.FaceAlpha = 0.5;
him.EdgeAlpha = 0.5;

maxval = 26; % somewhere around max of histogram
dotscale = 5;
scatter(good_pow_rest,dotscale+maxval+dotscale*rand(size(good_pow_diff)),15,'.','k')
scatter(good_pow_mvmt,maxval+dotscale*rand(size(good_pow_diff)),15,'.','r')
med_rest = median(good_pow_rest);
med_mvmt = median(good_pow_mvmt);
plot([med_rest med_rest],[dotscale+maxval 2*dotscale+maxval],'k')
plot([med_mvmt med_mvmt],[maxval dotscale+maxval],'r')