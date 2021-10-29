%%%%% Calculates the percentage of 0,1,2,3+ abx vs. met genes found per
%%%%% plasmid from prevalent and other STs and plots in stacked bar graph.
%%%%% Comment out line 15 to plot data from prevalent STs only; comment out 
%%%%% line 16 for all other STs.

close all, clear all, clc

% load in plasmid data filtered for transferrable plasmids only
Tplas = readtable("TableS3.xlsx");

% filter plasmid table for 4 prevalent STs
ST = str2double(Tplas.ST);

%%% comment out either line below (15 or 16) to plot respective data
% Tplas = Tplas(ismember(ST,[131 11 95 73]),:);
% Tplas = Tplas(~ismember(ST,[131 11 95 73]),:);

% calculate number of plasmids with no abx or met genes?
num_R_0 = sum(cellfun(@isempty,Tplas.resistance_drug_class));
num_M_0 = sum(Tplas.number_of_kegg_metabolism == 0);

% find all instances where resistance and metabolic genes exist, separately
Tplas2 = Tplas(~cellfun(@isempty,Tplas.resistance_drug_class),:);
Tplas3 = Tplas(Tplas.number_of_kegg_metabolism>0,:);

% determine number of resistance genes per plasmid
num_uniq_R = zeros(length(Tplas2.resistance_drug_class),1);
num_R = zeros(length(Tplas2.resistance_drug_class),1);
for q = 1:length(Tplas2.resistance_drug_class)
    if contains(Tplas2.resistance_drug_class{q},';')
        num_uniq_R(q) = length(unique(split(Tplas2.resistance_drug_class{q},';')));
        num_R(q) = length((split(Tplas2.resistance_drug_class{q},';')));
    else
        num_uniq_R(q) = 1;
        num_R(q) = 1;
    end
end

% determine number of metabolic genes per plasmid
ind=[];
kegg = readtable("TableS9.xlsx");
for q = 1:length(kegg.Column_Name)
    ind = [ind,find(strcmp(kegg.Column_Name{q},Tplas.Properties.VariableNames))];
end
Tkegg = Tplas3(:,ind);
num_uniq_M=[];
num_M=[];
for q = 1:height(Tplas3)
    num_uniq_M(q) = sum(table2array(Tkegg(q,:))>0);
    num_M(q) = sum(table2array(Tkegg(q,:)));
end

% concatenate gene counts per category (met or abx)
num_uniq_M = [num_uniq_M,zeros(1,num_M_0)];
num_uniq_R = [num_uniq_R;zeros(num_R_0,1)];

% determine number of 0,1,2, and 3+ and percentages for each gene type
% per plasmid 
tM = tabulate(num_uniq_M);
tR = tabulate(num_uniq_R);

str = 4;
tM(str,2:3) = sum(tM(str:end,2:3))
tM(str+1:end,:) = [];
tR(str,2:3) = sum(tR(str:end,2:3))
tR(str+1:end,:) = [];

%%%% generate stacked bar plot figure
figure; hold on
b1 = bar(1,[tR(:,3)]','stacked', 'linewidth',2.0); hold on
b2 = bar(2,[tM(:,3)]','stacked', 'linewidth',2.0)

set(gca,'xtick',[1 2], 'XTickLabel',{'Abx','Met'},...
    'fontsize',34,'linewidth',3.0,'XTickLabelRotation',90)
ylabel('% genes per plasmid')
b1(1).FaceColor = [229 247 255]./255;
b1(2).FaceColor = [102 212 255]./255;
b1(3).FaceColor = [0 160 229]./255;
b1(4).FaceColor = [0 110 153]./255;
b2(1).FaceColor = [255 229 229]./255;
b2(2).FaceColor = [255 0 0]./255;
b2(3).FaceColor = [204 0 0]./255;
b2(4).FaceColor = [127 0 0]./255;

% set legend for abx genes 
legend('0','1','2','3+','location','northoutside','orientation',...
    'horizontal')

set(gca,'tickdir','out')
set(gcf,'position',[560   425   360   483])
ylim([0 100]), set(gca,'ytick',[0 50 100])

%%%% likely a z test of proportions
% n=1775;
% p0 = .110986;
% sd = sqrt(p0 .* (1-p0) ./ n);
% z = (.196620 - .110986) ./ sd;

