%%%%% Calculates the frequency of abx vs. met genes per total plasmid genes
%%%%% and plots in log scale on histogram

close all, clear all, clc

% load in plasmid data filtered for transferrable plasmids only
Tplas = readtable("TableS3.xlsx");

% perform calculations for the log of gene types per total genes on plasmid
met = log(Tplas.number_of_kegg_metabolism./Tplas.num_genes);
abx = log(Tplas.number_of_resistance_genes./Tplas.num_genes);

%%%% generate and label histogram
figure; 
set(gcf,'position',[332   400   600   662])
p = panel();
p.pack(3,2);
p.de.margin = 2;
p.margin = [20 40 10 10];
p(1,1).select();
h1 = histogram(abx,80)
h1.FaceColor = [1, 178, 255]./255;
ylim([0 100])
set(gca,'FontSize',30,...
    'LineWidth',3.0,'TickDir','out','xtick',[])
line([mean(abx(abx>-inf)), mean(abx(abx>-inf))],[0 100],...
    'color','k','linestyle',':','linewidth',3.0)
box on
xlim([-6.2580   -1.6820])
p(3,1).select();
h2 = histogram(met,80)
h2.FaceColor = [255 25 63]./255 ;
ylim([0 100])
set(gca,'FontSize',30,...
    'LineWidth',3.0,'TickDir','out','ytick',[0 50 100],'xtick',[-6 -4 -2])
xlim([-6.2580   -1.6820])
line([mean(met(met>-inf)), mean(met(met>-inf))],[0 100],...
    'color','k','linestyle',':','linewidth',3.0)
box on
p.fontsize = 30;

% calculate percentage of MET and ABX per plasmid gene
pct_met_per_genes = (mean(Tplas.number_of_kegg_metabolism./...
    Tplas.num_genes))*100;
pct_abx_per_genes = (mean(Tplas.number_of_resistance_genes./...
    Tplas.num_genes))*100;

