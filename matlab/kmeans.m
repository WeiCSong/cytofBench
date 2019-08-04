clear;
% Kmeans for preprocessed data. When using, manually modify lines 3, 4, 13, 17, and 19
DataFile = {'colon'};
N = 20000;

R=[];
E=[];
for t=1:5
    %% random select N events
    tic
    
    % for preprocessed data
    fcsdat = readtable([pwd '/' DataFile{1} '/colon_20k/colon_20k_1.csv']);
    fcsdat = table2array(fcsdat(:, 2:end));
    [m,n]=size(fcsdat);
    
    rlabels=fcsdat(:,21); % select label column in X
    k=max(rlabels);
    X=fcsdat(:,1:19); % select marker columns in X
    X=asinh(X./5);
    %% kmeans clustering
    opts = statset('Display','final','MaxIter',1000);
    [Idx,C,sumD,D]=kmeans(X,k,'Options',opts);
    subpop{t}.X=X;
    subpop{t}.Idx=Idx;%label of each selected event
    subpop{t}.C=C;%k*n matrix,K个类的质心位置
    subpop{t}.sumD=sumD;%
    subpop{t}.D=D;
    subpop{t}.rlabels=rlabels;
    %% Evaluation computation
    res=bestMap(rlabels,Idx);
    ac(t) = length(find(rlabels == res))/length(rlabels);
    nmi(t)=compute_NMI(rlabels,Idx);
    P=rlabels';
    Q=Idx';
    FMeasure(t) = Fmeasure(P,Q);
    
    %% distribution of cluster size
    N1=length(rlabels);
    rl=cell(k,1);
    el=cell(k,1);
    for i=1:k
        rl{i}=find(rlabels==i);
        rs(i)=length(rl{i})/N1;%真实的每个类大小比例
        el{i}=find(res==i);
        es(i)=length(el{i})/N1;%算法的每个类大小比例
    end
    R=[R;rs];
    E=[E;es];
    Time{t} = toc;
    
    mkdir([pwd '/' DataFile{1}])
    mkdir([pwd '/' DataFile{1} '/' num2str(N)])
    save( [pwd '/' DataFile{1} '/' num2str(N) '/' DataFile{1} '.mat'],'subpop')
end
%% compute average measurement
Mac=mean(ac)
Sac=std(ac)
Mnmi=mean(nmi)
Snmi=std(nmi)
MFMeasure=mean(FMeasure)
SFMeasure=std(FMeasure)

%% plot cluster size distribution
MR=mean(R);%10次均值
SR=std(R);%10次标准差
ME=mean(E);%10次均值
SE=std(E);%10次标准差
clusterSize=[MR' ME'];
sstd=[SR' SE'];
bar(clusterSize);
hold on;
errorbar(clusterSize,sstd,'LineStyle','none')
title(DataFile{1},'FontSize',12)
xlabel('# cluster number','FontSize',12);
ylabel('# cluster size','FontSize',12);
legend('Gating','kmeans')
saveas(gcf, [pwd '/' DataFile{1} '/' num2str(N) '/' DataFile{1} '_ClusterSize.png']);

true=MR';
pred=ME';
figure,scatter(true,pred,50,'filled')
box on, grid on
xlabel('True frequency'),ylabel('Predicted frequency')
title(DataFile{1})
lsline
saveas(gcf, [pwd '/' DataFile{1} '/' num2str(N) '/' DataFile{1} '_Frequency.png']);

for t=1:5
    disp(['Runtime: ' num2str(Time{t}) ' seconds'])
end
