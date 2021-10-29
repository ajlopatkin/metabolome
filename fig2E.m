%%%%% Calculates the percentage of genes from kegg metabolic categories 
%%%%% and known resistance drug classes and plots on a horizontal bar
%%%%% graph

close all, clear all, clc

% load in plasmid data filtered for transferrable plasmids only
Tplas = readtable("TableS3.xlsx");

% load in met categories and calculate percentage
mets = readtable("TableS9.xlsx");met = mets.Column_Name;
Tmet = [];
for q = 1:length(met)
    Tmet(:,end+1) = double(Tplas.(met{q})>0);
end
totalM = [sum(Tmet)./(sum(sum(Tmet))).*100];

% load in drug categories and calculate percentage
drugs = readtable("TableS8.xlsx");  
drugs([3,6,7,10,14,15],:)=[]; % filter out categories not of interest
drug = drugs.drug_class;
all_res_genes = split(join(Tplas.resistance_drug_class,';'),';');
num_R = zeros(length(drug),1);
for q = 1:length(drug)
    num_R(q) = sum(strcmp(all_res_genes,drug{q}));
end
totalR = num_R./sum(num_R).*100;

%%%% generate horizontal bar graphs on a panel
figure; 
set(gcf,'position',[332   193   854   662])
p = panel();
p.pack(2,2);
p.de.margin = 2;
p.margin = [30 50 10 10];
h1=p(1,2).select();

barh(totalM,'linewidth',2.0,'facecolor',[255 25 63]./255)
xlim([0 30])
set(gca,'ytick',[1:length(totalM)],...
    'YTickLabel',mets.Legend_Name,...
    'YTickLabelRotation',0,...
    'FontSize',30,...
    'LineWidth',3.0,...
    'xtick',[],'TickDir','out')
h1=p(2,2).select();

barh(totalR,'linewidth',2.0,'facecolor',[1, 178, 255]./255)
xlim([0 30])
set(gca,'ytick',[1:length(totalR)],...
    'YTickLabel',drugs.legend_names,...
    'YTickLabelRotation',0,...
    'XTick',[0 10 20 30],...
    'FontSize',30,...
    'LineWidth',3.0,'TickDir','out')
p.fontsize = 34; 
