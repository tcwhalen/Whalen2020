function plotFracs(fracs,cols,shapes,titletext,xlabels)
% Tim C Whalen, last edited July 2020
% Helper for plotting oscillation fractions in figures 3, 7 and 9 from
% Whalen et al. 2020 JNeurophys

cols_dark = cols./1.3;
nconds = size(fracs,1);
for c = 1:nconds
    fr = fracs{c};
    np = length(fr);
    me = mean(fr);
    plot([-.1 .1]+c,[me me],'Color',1-(1-cols(c,:))/1.5,'LineWidth',2)
    scatter(c+((0:.5:(np-1)/2)-np/4)/10, fr, 15,cols(c,:),shapes(c),'MarkerEdgeColor',cols_dark(c,:),'MarkerFaceColor',cols(c,:))
end
ylim([0 1])
xlim([.5 4.5])
set(gca, 'xtick',1:nconds,'XTickLabels',xlabels,'ytick', [0 1]) %,'xcolor','w')
set(gca, 'xtick',[],'ytick', [0 1])
title(titletext)
makeNice(gca)
