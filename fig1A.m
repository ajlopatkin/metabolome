%%%% Calculates the percentage of each of the four prevalent STs we
%%%% highlight in our dataset and other STs and plots on bar graph

close all, clear all, clc

% load in plasmid genome table and filter out those with undefined STs
T_plasmid = readtable("TableS10.xlsx");
T_plasmid(find(strcmp(T_plasmid.ST,'-')),:) = [];

% load in table with data collected from literature search
T_litsearch = readtable("TableS1.xlsx");

% define STs in dataset
dataSTs = T_plasmid.ST;
litSTs = T_litsearch.LITEARCH_ST_label;
C = setdiff(str2double(dataSTs),litSTs);
ind = [];
for q = 1:length(T_plasmid.ST)
    currST = str2double(T_plasmid.ST{q});
    if isempty(find(C == currST))
        ind(end+1) = q;
    end
end
T_plasmid = T_plasmid(ind,:);

% determine prevalent STs vs. other STs frequency in dataset 
prevalent_ST = T_litsearch.LITEARCH_ST_label(T_litsearch.LITEARCH_ST_pct>0.5);
tab = tabulate(T_plasmid.ST);
STs = str2double(tab(:,1));
freq = [tab{:,3}]';
tot = [];
st_labs = [];
for q = 1:length(prevalent_ST)
    curSTS = freq(STs == prevalent_ST(q));
    if isempty(curSTS)
        curSTS = 0;
    end
    
    if curSTS > 7
    tot = [tot,curSTS];
    st_labs = [st_labs,prevalent_ST(q)];
    end
end
[tot i] = sort(tot,'descend')
st_labs = st_labs(i);

other = [zeros(1,length(tot)-1),100-sum(tot)];
tot = [tot;other];

%%%% generate figure
figure; hold on
b1 = barh([1,2],tot,'stacked','linewidth',2.0)
b2 = barh([2],other,'stacked','facecolor',[253 253 150]./255,'linewidth',2.0)
b1(1).FaceColor = [239 190 125]./255;
b1(2).FaceColor = [202 231 193]./255;
b1(3).FaceColor = [201 214 232]./255;
b1(4).FaceColor = [255 215 250]./255;

xlim([0 100])
legend([b1(1) b1(2) b1(3) b1(4)],...
    {string(st_labs(1)),string(st_labs(2)),string(st_labs(3)),...
    string(st_labs(4))},'Location','northoutside','orientation','horizontal','fontsize',26,'linewidth',2.0)
set(gca,'ytick',[1 2],'yticklabel',{"Prevalent","Other"},...
    'yTickLabelRotation',0,'fontsize',40,'linewidth',4.0)
