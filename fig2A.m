%%%%% Calculates the percentage of metabolic genes, resistance genes, and 
%%%%% all other genes and plots on pie chart
close all, clear all, clc

% load in plasmid data filtered for transferrable plasmids only
Tplas = readtable("TableS3.xlsx");

% perform gene count and percentage calculations
num_M = sum(Tplas.number_of_kegg_metabolism);
num_R = sum(Tplas.number_of_resistance_genes);
num_genes = sum(Tplas.num_genes);
pctR = [num_R./num_genes*100];
pctM = [num_M./num_genes*100];
pctO = 100 - (pctM+pctR);

%%%% generate and label pie chart
figure; 
h=pie([pctR,pctM,pctO])
h(1).FaceColor = [1, 178, 255]./255;
h(2).FontSize = 30;
h(2).Position = [.03 1.0958 0];
h(2).Color = [1, 178, 255]./255;
h(3).FaceColor = [255 25 63]./255 ;
h(4).Color = [255 25 63]./255
h(4).FontSize = 30;
h(4).Position = [-.48 1 0]
h(5).FaceColor = [.6 .6 .6]
h(6).FontSize = 30;
h(6).Color = [.6 .6 .6];
legend('Abx','Met','Other','fontsize',30,'location','eastoutside','orientation','vertical')
set(gca,'linewidth',3.0)

h(2).String = string(round(pctR,1)) + "%";
h(4).String = string(round(pctM,1)) + "%";
h(6).String = string(round(pctO,1)) + "%";