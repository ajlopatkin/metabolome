%%%%% Calculates the percentage of metabolic L2 and resistance type genes
%%%%% and plots on bar graph

close all, clear all, clc
addpath('freezeColors/')
addpath('cbrewer/')

% load in plasmid data and filter for transferrability
Tplas = readtable("TableS3.xlsx");

% perform met gene vs. abx gene per total gene calculations
met = (Tplas.number_of_kegg_metabolism./Tplas.num_genes)*100;
abx = (Tplas.number_of_resistance_genes./Tplas.num_genes)*100;
oth = ((Tplas.num_genes - ...
    (Tplas.number_of_kegg_metabolism+Tplas.number_of_resistance_genes))./Tplas.num_genes)*100;
[h i] = sort(met,'ascend');
mat = [met(i),abx(i),oth(i)];

% define inc groups
Tinc = readtable("TableS12.xlsx"); incs=Tinc.INC; % replace incs with: incs = unique(strsplit(strjoin(Tplas.PlasFinder_SIMPLIFIED, ';'), ';'))'
mat_inc = zeros(height(Tplas),length(incs));
for q = 1:length(incs)
    mat_inc(:,q) = contains(Tplas.PlasFinder_SIMPLIFIED,incs{q});
end


%%%% generate figure
p = panel();
p.pack('h', {1/3 0.5/3 1.5/3})
p.de.margin = 8;
p.margin = [10 30 2 10];
set(gcf,'position',[332   400   600   662])

% (Panel 1) ST heatmap
p(1).pack('h',{1/8 1/8 []});
p(1,1).select()
ST = str2double(Tplas.ST(i));
ST(~(ST == 131 | ST == 73 | ST == 95 | ST == 11)) = 0;
ST(ST == 131) = 1;
ST(ST == 73) = 2;
ST(ST == 95) = 3;
ST(ST == 11) = 4;
[cmap] = flipud(cbrewer('qual', 'Pastel2', 6)); cmap(5,:) = [];
imagesc(ST)
ylim([1 height(Tplas)])
set(gca,'xtick',[],'ytick',[],'linewidth',2.0), box on
colormap(cmap), freezeColors

% (Panel 2) mobility type heatmap 
p(1,2).select()
MOB = (Tplas.PredictedMobility(i));
MOB_new = zeros(height(Tplas),1);
MOB_new(strcmp(MOB,'Mobilizable')) = 1;
[cmap] = flipud(cbrewer('qual', 'Dark1', 2)); 
imagesc(MOB_new)
ylim([1 height(Tplas)])
set(gca,'xtick',[],'ytick',[],'linewidth',2.0), box on
colormap([169 169 169; 43 45 47]./255), freezeColors

% (Panel 3) stacked bar graph
p(1,3).select()
b = barh(mat,'stacked')
b(1).FaceColor = [255 25 63]./255; 
b(2).FaceColor = [1, 178, 255]./255;
b(3).FaceColor = [0 0 0];
set(gca,'xscale','log','fontsize',12,'ytick',[],...
    'linewidth',2.0,'TickDir','out','TickLength',[0 0]), box on
ylim([1 height(Tplas)]), xlim([10^-2 10^2])

% (Panel 4) bar graph
p(2).select()
barh(Tplas.num_genes(i),'FaceColor',[.6 .6 .6])
set(gca,'fontsize',12,'ytick',[],...
    'linewidth',2.0,'xtick',[0 250 500],'TickDir','out','TickLength',[0 0]), box on
ylim([1 height(Tplas)])

% (Panel 5) inc heatmap
p(3).select()
imagesc(mat_inc)
ylim([1 height(Tplas)]), xlim([0.5 length(incs)+0.5])
set(gca,'xtick',1:size(mat_inc,2),'xticklabel',incs,'XTickLabelRotation',90,...
    'fontsize',14,'ytick',[],'linewidth',2.0,'TickDir','out','TickLength',[0 0]), box on
colormap([1 1 1;0 0 0])
