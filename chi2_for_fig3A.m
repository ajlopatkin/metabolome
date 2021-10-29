
%%%%% Individual genes from significant metabolic and drug classes
%%%%% are compared against one another using the fisher exact test.
%%%%% Tables are generated demonstrating significant relationships 
%%%%% between each gene.


clear, close all

% load in plasmid data filtered for transferrable plasmids only
Tplasmid = readtable("TableS3");
Tgenedata = readtable("TableS12.xlsx");

% generate list of unique MET and ABX genes
allgenes = split(join(cellstr((Tplasmid.all_kegg_metabolism_genes)),';'),';');
unique_genes_MET = unique((allgenes(~cellfun(@isempty,allgenes))));
allgenes = split(join(cellstr((Tplasmid.resistant_genes)),';'),';');
unique_genes_ABX = unique((allgenes(~cellfun(@isempty,allgenes))));
genes = [unique_genes_ABX;unique_genes_MET];
type = [repelem("ABX",length(unique_genes_ABX),1);repelem("MET",length(unique_genes_MET),1)];
unique_drug_class = rmmissing(unique(strsplit(strjoin(Tgenedata.RESISTANCE, ';'), ';')))';
class = cell(length(genes), 1);

% determine KEGG metabolism class per MET gene
for g = 1: length(genes)
    
    row = find(strcmp(Tgenedata.gene, genes(g)));
    if Tgenedata.is_kegg_carb(row) == 1
        class(g) = {'CARB'};
    end
    
    if Tgenedata.is_kegg_energy(row) == 1
        class(g) = {'ENERGY'};
    end
    
    if Tgenedata.is_kegg_lipid(row) == 1
        class(g) = {'LIPID'};
    end
    
    if Tgenedata.is_kegg_nucleotide(row) == 1
        class(g) = {'NUCLEOTIDE'};
    end
    
    
    if Tgenedata.is_kegg_aminoacid(row) == 1
        class(g) = {'AMINO ACIDS'};
    end
    
    
    if Tgenedata.is_kegg_other_aminoacids(row) == 1
        class(g) = {'OTHER AMINO ACIDS'};
    end
    
    if Tgenedata.is_kegg_glycan_biosynthesis(row) == 1
        class(g) = {'GLYCAN BIOSYNTHESIS'};
    end
    
    
    if Tgenedata.is_kegg_cofactors(row) == 1
        class(g) = {'COFACTORS'};
    end
    
    
    if Tgenedata.is_kegg_terpenoids(row) == 1
        class(g) = {'TERPENOIDS'};
    end
    
    
    if Tgenedata.is_kegg_sec(row) == 1
        class(g) = {'SECONDARY'};
    end
    
    if Tgenedata.is_kegg_xenobiotics(row) == 1
        class(g) = {'XENOBIOTICS'};
    end
    
    if Tgenedata.is_kegg_not_included(row) == 1
        class(g) = {'NOT INCLUDED'};
    end
    
    
end

% determine drug class per ABX gene
row2 =  find(cellfun(@isempty, class));
all_res_genes = genes(row2);

for gg = 1:length(all_res_genes)
    
    row3 = find(strcmp(Tgenedata.gene, all_res_genes(gg)));
    drug_class = Tgenedata.RESISTANCE(row3(1));
    
    row4 = find(strcmp(genes, all_res_genes(gg)));
    class(row4) = drug_class;
    
end

% add information to unique MET and ABX gene list
if ~isfile("TableS16.xlsx")
    writetable(array2table([type,genes,class],'VariableNames',{'type','genes', 'class',;}),"TableS14_unique_gene_list.xlsx")
end

% remove drug and metabolism classes that are not of interest
Tgene = readtable("TableS16.xlsx");
Tgene(strcmp(Tgene.class,''),:) = [];
Tgene(strcmp(Tgene.class,'BLEOMYCIN'),:) = [];
Tgene(strcmp(Tgene.class,'FOSFOMYCIN'),:) = [];
Tgene(strcmp(Tgene.class,'OTHER AMINO ACIDS'),:) = [];
Tgene(strcmp(Tgene.class,'PEPTIDE'),:) = [];
Tgene(strcmp(Tgene.class,'RIFAMYCIN'),:) = [];

% define updated version of gene and gene type list
genes = Tgene.genes;
type = Tgene.type;


% define variables to search
METABOLISM_GENES = Tplasmid.all_kegg_metabolism_genes;
RESISTANT_GENES = Tplasmid.resistant_genes;

% make new table to fill in with p values, this is will be a square
% matrix where the diagonal compares the gene against itself
T_chisq = [table('size',[length(genes),length(genes)],'VariableTypes',...
    repelem({'double'},1,length(genes)),'VariableNames',genes)];

% initiate collection matrix
gene1 = [];
gene2 = [];
gene1_type = [];
gene2_type = [];
pval = [];
conting = [];
odds_rat = [];

% loop through every gene twice
for x = 1:length(genes)
    
    if ~mod(x,50)
        disp("Loop " + x + " out of " + length(genes))
    end
    
    for y = 1:length(genes)
        
        % determine plasmids with gene 1
        cur_gene_1 = genes{x};
        type_1 = type{x};
        if strcmp(type_1,'MET')
            ind_1 = contains(METABOLISM_GENES,cur_gene_1);
        else
            ind_1 = contains(RESISTANT_GENES,cur_gene_1);
        end
        
        % determine plasmids with gene 2
        cur_gene_2 = genes{y};
        type_2 = type{y};
        if strcmp(type_2,'MET')
            ind_2 = contains(METABOLISM_GENES,cur_gene_2);
        else
            ind_2 = contains(RESISTANT_GENES,cur_gene_2);
        end
        
        % when the two genes are the same, break and go to next loop
        if strcmp(cur_gene_2,cur_gene_1)
            T_chisq(x,y:end) = array2table(0);
            break
        end
        
        
        % calculate chisquare
        [tbl c p] = crosstab(ind_1,ind_2);
        
        [h,p,stats] = fishertest(tbl);
        
        T_chisq(x,y) = array2table(p);
        
        
        % calculate chisquare
        [tbl c p] = crosstab(ind_1,ind_2);
        [h,p,stats] = fishertest(tbl);
        T_chisq(x,y) = array2table(p);
        
        % collect names and types of significant pairs of genes
        if (p < .05)
            gene1 = [gene1;string(cur_gene_1)];
            gene1_type = [gene1_type;string(type_1)];
            gene2 = [gene2;string(cur_gene_2)];
            gene2_type = [gene2_type;string(type_2)];
            pval = [pval;p];
            conting = [conting;reshape(tbl,1,4)];
            odds_rat = [odds_rat;stats.OddsRatio];
        end
        
    end
end

% add columns to gene1 table for gene comparisons by for loop below
has_gene1 = zeros(height(gene1),1);
has_gene2 = zeros(height(gene1),1);
gene_both = zeros(height(gene1),1);
gene1only = zeros(height(gene1),1);
gene2only = zeros(height(gene1),1);
gene_neither = zeros(height(gene1),1);


for z = 1:length(gene1)
    
    gene_columns = [Tplasmid.all_kegg_metabolism_genes, Tplasmid.resistant_genes];
    check1 = [gene1(z)];
    check2 = [gene2(z)];
    
    gene1_check = contains(gene_columns(:,1), check1) + contains(gene_columns(:,2), check1);
    gene2_check = contains(gene_columns(:,1), check2) + contains(gene_columns(:,2), check2);
    
    has_gene1(z) = sum(gene1_check);
    has_gene2(z) = sum(gene2_check);
    
    gene_both_row = zeros(height(gene1_check),1);
    gene1only_row = zeros(height(gene1_check),1);
    gene2only_row = zeros(height(gene1_check),1);
    gene_neither_row = zeros(height(gene1_check),1);
    
    for zz = 1: length(gene1_check)
        % determine if plasmid has both
        gene_both_row(zz) = ((sum(gene1_check(zz) + gene2_check(zz)) == 2));
    
        % determine if plasmid has only gene 1
        gene1only_row(zz) = (gene1_check(zz) == 1 && gene2_check(zz) == 0);
      
        
        % determine if plasmid has only gene 2
        gene2only_row(zz) = (gene1_check(zz) == 0 && gene2_check(zz) == 1);
      
        
        % determine if plasmid has neither gene 1 or gene 2
        gene_neither_row(zz) = (gene1_check(zz) == 0 && gene2_check(zz) == 0);
        
        
    end
     gene_both(z) = sum(gene_both_row);
     gene1only(z) = sum(gene1only_row);
     gene2only(z) = sum(gene2only_row);
     gene_neither(z) = sum(gene_neither_row);
    
end

% save data
T = [array2table(type),array2table(genes)];
T_signif_genes = table(gene1,gene1_type,gene2,gene2_type,pval,odds_rat,has_gene1,has_gene2, gene_both, gene1only, gene2only,gene_neither);

writetable(T_signif_genes,"TableS14.xlsx")
T_chisq_save = [T,T_chisq];
writetable(T_chisq_save,"TableS13.xlsx")

%%%% generate raw chi square heatmap
p_square = table2array(T_chisq(:,3:end));
figure; imagesc(p_square)
caxis([0 1])
set(gca,'ytick',[1:height(T_chisq)],'YTickLabel',type)