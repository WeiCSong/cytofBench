%% Read the Data and Preprocess
clear
FolderName = {'colon'};
N = 20000;

for t = 1:5
DataTable = readtable([pwd '/colon/colon_' num2str(N/1000) 'k/colon_' num2str(N/1000) 'k_' num2str(t) '.csv']);

% remove unneeded columns
%{
DataTable.Time=[];
DataTable.Cell_length=[];
DataTable.DNA1=[];
DataTable.DNA2=[];
DataTable.Viability=[];
DataTable.file_number=[];
DataTable.event_number=[];
DataTable.subject=[];
%}
DataTable.Var1=[];

% Separate Data points and Labels
Labels = DataTable.label;
DataTable.label=[];
Data = table2array(DataTable);

clear DataTable

% clear NotDebrisSinglets
Data(strcmp('NotDebrisSinglets',Labels),:)=[];
Labels(strcmp('NotDebrisSinglets',Labels))=[];



%% run LDA Classifier with percent cross-validation

CVO = cvpartition(Labels,'k',5);
Accuracy = zeros(CVO.NumTestSets,1);
nmi=zeros(CVO.NumTestSets,1);
FMeasure=zeros(CVO.NumTestSets,1);
training_time = zeros(CVO.NumTestSets,1);
testing_time = zeros(CVO.NumTestSets,1);
%CellTypes = {'','','',};
numCellTypes = unique(Labels);
ConfusionMat = zeros(length(numCellTypes));
for i = 1:CVO.NumTestSets
    trIdx = CVO.training(i);
    teIdx = CVO.test(i);
    tic
    classificationLDA = fitcdiscr(...
        Data(trIdx,:), ...
        Labels(trIdx));
    training_time(i)=toc;          %in seconds
    
    tic
    Predictor = predict(classificationLDA,Data(teIdx,:));
    %Predictor = string(Predictor)
    testing_time(i)=toc;           %in seconds
    rlabels=Labels(teIdx);
    Result{i}.rlabels=rlabels;
    Result{i}.labels=Predictor;
    Accuracy(i) = length(find(rlabels == Predictor))/length(rlabels);
    nmi(i)=compute_NMI(rlabels,Predictor);
    P=rlabels';
    C=Predictor';
    %testtest = Fmeasure(P,C);
    FMeasure(i) = Fmeasure(P,C);
    ConfusionMat = ConfusionMat + confusionmat(rlabels,Predictor,'order',numCellTypes);
    trdataID=find(trIdx);
    Result{i}.trdataID=trdataID;
    tedataID=find(teIdx);
    Result{i}.tedataID=tedataID;
    Result{i}.teData=Data(teIdx,:);
    Result{i}.trData=Data(trIdx,:);
end
save([pwd '/' FolderName{1} '/' num2str(N) '/ldaResult_' FolderName{1} '_' num2str(t) '.mat'],'Result')
Total_time(t) = sum(training_time)+sum(testing_time);
training_time(t) = mean(training_time);
testing_time(t) = mean(testing_time);
cvAcc(t) = mean(Accuracy)*100;
cvSTD(t) = std(Accuracy)*100;
cvnmi(t) = mean(nmi);
cvstdnmi(t) = std(nmi);
cvfm(t) = mean(FMeasure);
cvstdfm(t) = std(FMeasure);
disp(['LDA Accuracy = ' num2str(cvAcc(t)) ' ' char(177) ' ' num2str(cvSTD(t)) ' %'])
disp(['LDA NMI = ' num2str(cvnmi(t)) ' ' char(177) ' ' num2str(cvstdnmi(t)) ' %'])
disp(['LDA FMeasure = ' num2str(cvfm(t)) ' ' char(177) ' ' num2str(cvstdnmi(t)) ' %'])
clear i Predictor classificationLDA trIdx teIdx CVO 

%% Performance evaluation

% F1 measure
Precision = diag(ConfusionMat)./sum(ConfusionMat,1)';
Recall = diag(ConfusionMat)./sum(ConfusionMat,2);
F1measure = 2 * (Precision.*Recall)./(Precision+Recall);
MedianFmeasure = median(F1measure);
Subset_size = sum(ConfusionMat,2);
WeightedFmeasure = (Subset_size./size(Data,1))'*F1measure;

disp(['Median F1-score = ' num2str(MedianFmeasure)])
figure,scatter(log10(Subset_size),F1measure,100,'filled'),title(FolderName{1})
xlabel('Log10(population size)'),ylabel('F1-score'),box on, grid on
saveas(gcf, [pwd '/' FolderName{1} '/' num2str(N) '/Performance_Evaluation_' num2str(t) '.png']);

%% Population Frequency

True_Freq = sum(ConfusionMat,2)./sum(sum(ConfusionMat));
Predicted_Freq = sum(ConfusionMat,1)'./sum(sum(ConfusionMat));
Max_Freq_diff = max(abs(True_Freq-Predicted_Freq))*100;

disp(['delta_f = ' num2str(Max_Freq_diff)])
figure,bar([True_Freq*100 Predicted_Freq*100])
%xticklabels(CellTypes)
%xtickangle(90)
set(gca,'FontSize',15)
legend({'True','Predicted'},'FontSize',15)
legend show
ylabel('Freq. %'),title(FolderName{1})
saveas(gcf, [pwd '/' FolderName{1} '/' num2str(N) '/Population_Frequency_' num2str(t) '.png'])

%% Population Frequency scatter plot
%{
X=log(True_Freq*100);
Y=log(Predicted_Freq*100);
figure,scatter(X,Y,50,'filled')
box on, grid on
xlabel('Log(True frequency %)'),ylabel('Log(Predicted frequency %)')
title('AML')
for k=1:length(CellTypes)
    text(X(k),Y(k),CellTypes{k})
end
lsline
text(0,0,['R = ' num2str(corr(X,Y))])
saveas(gcf, [pwd '/' FolderName{1} '/' num2str(N) '/Population_Frequency_scatter_log_' num2str(t) '.png']);
%}

X=True_Freq;
Y=Predicted_Freq;
figure,scatter(X,Y,50,'filled')
box on, grid on
xlabel('True frequency'),ylabel('Predicted frequency')
title(FolderName{1})
lsline
text(0,0,['R = ' num2str(corr(X,Y))])
saveas(gcf, [pwd '/' FolderName{1} '/' num2str(N) '/Population_Frequency_scatter_' num2str(t) '.png']);
disp(['Training time: ' num2str(training_time(t))])
disp(['Testing time: ' num2str(testing_time(t))])
disp(['Total time: ' num2str(Total_time(t))])
end
Mac=mean(cvAcc)
Sac=std(cvAcc)
Mnmi=mean(cvnmi)
Snmi=std(cvnmi)
MFMeasure=mean(cvfm)
SFMeasure=std(cvfm)
save([pwd '/' FolderName{1} '/' num2str(N) '/all_results.mat'])
