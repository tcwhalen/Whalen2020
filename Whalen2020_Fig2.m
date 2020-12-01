% Tim C Whalen, last edited July 2020
% Makes Whalen 2020 JNeurophys Figure 2B-C, examples of true and false
% positive oscillation detections
load('Whalen2020_data_fig2','data');
data.osc = struct();
data.osc.toPlot=1;

% examples to plot - units 1, 2, 1 and 1 from files 1-4.
data.osc = struct('examples',[1 2; 2 2; 3 1; 4 1]);

[ data ] = renewalPSD_phaseShift_batch(data);
% Note that this function only displays red dots for oscillations that
% pass both criteria. To see false positives (i.e. passed PSD but not phase 
% shift) for examples 3 and 4, compare data.osc.has_osc (1 if unit has
% oscillation based on PSD and phase shift) to data.has_osc_old (1 if unit
% has oscillation based solely on PSD)