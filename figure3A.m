%%%%% Co-occuring abx genes are tabulated for each significant metabolic 
%%%%% gene and these associations are represented on a panel plot with 
%%%%% indications of metabolic class and number of significantly 
%%%%% disassociated genes (in black) per pair represented.


close all, clear all, clc

% load in plasmid-level dataset, gene-level dataset, and chi square data

Tplas = readtable("TableS3.xlsx");
Tgenes = readtable("TableS12.xlsx");
T = readtable("TableS14.xlsx"); 

% all metabolic genes the co-occur with resistance genes
types = strcat(T.gene1_type,T.gene2_type);
T_both = T((contains(types,'MET')&contains(types,'ABX')),:);
tab = tabulate(T_both.gene1);
tab(strcmp(tab(:,1),'fab'),:) = [];% remove due to redundancy(fabG present)
[~,i] = sort(cell2mat(tab(:,2)),'descend');
tabM = tab(i,:);
tabM(:,4) = cell(size(tabM,1),1);
Tkegg = readtable("TableS16.xlsx"); %currently saved in sheet 2 of Table S9
kegg_IDs = Tkegg.KEGG_ID;


% collect the metabolism categories for each gene
mat = zeros(size(tabM,1),length(Tkegg.KEGG_NAME));
count = 1;
for qq = 1:length(Tkegg.KEGG_ID)
    for q = 1:size(tabM,1)
        Tcurr = Tgenes(find(contains(Tgenes.gene,tabM{q,1})),:);
        tabM{q,4} = Tcurr.kegg_metabolism_class{1};
        curr_kegg = contains(Tcurr.kegg_metabolism_class{1},num2str(Tkegg.KEGG_ID(qq)));
        if curr_kegg
            mat(q,qq) = count;
        end
    end
    count = count+1;
end

% filter matrix by removing met genes with less than 2 co-occurring abx genes
ind = sum(mat,2)<1;
mat(ind,:) = [];
tabM(ind,:) = [];


% determine the number of significantly disassociated genes per gene (based
% on odds ratios (<1))
top_gene1 = tabM(:, 1);

for q = 1: length(top_gene1)
    row = find(strcmp(T_both.gene1, top_gene1(q)));
    oddsrat = T_both.odds_rat(row);
    num_disassoc = zeros(length(oddsrat), 1);
    for qq = 1: length(oddsrat)
        num_disassoc(qq) = (oddsrat(qq)<1);
    end
    tabM(q, 5) = {sum(num_disassoc)};
end


% flip matrix to orient bar plot and heatmap in same direction
mat = flipud(mat); 

% generate figure
p = panel();
p.pack('h', {1/3 2/3})
p.de.margin = 2;
p.margin = [30 60 2 10];
set(gcf,'position',[200   300   500   800])

% bar plot
p(1).select()
barh(flipud(cell2mat(tabM(:,2))),'facecolor',[.6 .6 .6],'linewidth',2.0); hold on
barh(flipud(cell2mat(tabM(:,5))),'facecolor',[0 0 0],'linewidth',2.0)
ylim([.5 size(tabM,1)+.5])
set(gca,'fontsize',30,'tickdir','out','linewidth',2.0,...
    'ytick',[1: size(tabM,1)],'yticklabel',flipud(tabM(:,1))), ylim([0.5 size(tabM,1)+.5])

% heatmap
p(2).select()
imagesc((mat))
xlim([.5 6.5]), ylim([.5 size(tabM,1)+.5])
set(gca,'ytick',[],'xtick',[1:6],...
    'XTickLabel',Tkegg.KEGG_NAME,'XTickLabelRotation',90,...
    'fontsize',40,'tickdir','out','linewidth',2.0), box on
cmap = [1 1 1;cbrewer('qual', 'Dark2', 7)];
cmap(end-2,:) = [];
colormap(cmap), hold on
for q = 1:size(tabM,1)+1
    plot([q-.5 q-.5],[0 size(tabM,1)+1],'k','linewidth',2.0)
    plot([0 7],[q-.5 q-.5],'k','linewidth',2.0)
end

% save table with number of associations vs. dissociations per gene
T_genes = table();
T_genes.gene_name = tabM(:,1);
T_genes.num_assoc = cell2mat(tabM(:,2))-cell2mat(tabM(:,5));
T_genes.num_dissoc = tabM(:,5);
writetable(T_genes,'TableS15_sheet3.xlsx');