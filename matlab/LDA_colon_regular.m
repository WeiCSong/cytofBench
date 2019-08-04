%% Read the Data and Preprocess
clear
FolderName = {'colon'};
mkdir([pwd '/' FolderName{1} '/Regular/Results'])

DataTable = readtable([pwd '/colon/colon.csv']);

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
Samples = DataTable.sample;
DataTable.label=[];
DataTable.sample=[];
Data = table2array(DataTable);

Samples_train = logical(zeros(size(Samples)));
Samples_train(Samples >= 17 & Samples <= 24) = 1;
Samples_test = logical(zeros(size(Samples)));
Samples_test(Samples < 17 | Samples > 24, :) = 1;


clear DataTable

% clear NotDebrisSinglets
Data(strcmp('NotDebrisSinglets',Labels),:)=[];
Labels(strcmp('NotDebrisSinglets',Labels))=[];

% Apply arcsinh5 transformation
Data=asinh((Data-1)/5);

%% run LDA Classifier with 5-fold cross-validation

% CVO = cvpartition(Labels,'k',2);
Accuracy = zeros(1,1);
nmi=zeros(1,1);
FMeasure=zeros(1,1);
training_time = zeros(1,1);
testing_time = zeros(1,1);
%CellTypes = {'','','',};
numCellTypes = unique(Labels);
ConfusionMat = zeros(length(numCellTypes));
for i = 1:1
    trIdx = Samples_train;
    teIdx = Samples_test;
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
    FMeasure(i) = Fmeasure(P,C);
    ConfusionMat = ConfusionMat + confusionmat(rlabels,Predictor,'order',numCellTypes);
    trdataID=find(trIdx);
    Result{i}.trdataID=trdataID;
    tedataID=find(teIdx);
    Result{i}.tedataID=tedataID;
    Result{i}.teData=Data(teIdx,:);
    Result{i}.trData=Data(trIdx,:);
end
save([pwd '/' FolderName{1} '/Regular/Results/ldaResult_colon_regular.mat'],'Result')
Total_time = sum(training_time)+sum(testing_time);
training_time = mean(training_time);
testing_time = mean(testing_time);
cvAcc = mean(Accuracy)*100;
cvSTD = std(Accuracy)*100;
cvnmi = mean(nmi);
cvstdnmi = std(nmi);
cvfm = mean(FMeasure);
cvstdfm = std(FMeasure);
disp(['LDA Accuracy = ' num2str(cvAcc) ' ' char(177) ' ' num2str(cvSTD) ' %'])
disp(['LDA NMI = ' num2str(cvnmi) ' ' char(177) ' ' num2str(cvstdnmi) ' %'])
disp(['LDA FMeasure = ' num2str(cvfm) ' ' char(177) ' ' num2str(cvstdnmi) ' %'])
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
saveas(gcf, [pwd '/' FolderName{1} '/Regular/Results/Performance_Evaluation.png']);

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
saveas(gcf, [pwd '/' FolderName{1} '/Regular/Results/Population_Frequency.png'])

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
saveas(gcf, [pwd '/' FolderName{1} '/Regular/Results/Population_Frequency_scatter.png']);
disp(['Training time: ' num2str(training_time)])
disp(['Testing time: ' num2str(testing_time)])
disp(['Total time: ' num2str(Total_time)])

Mac=mean(cvAcc)
Sac=std(cvAcc)
Mnmi=mean(cvnmi)
Snmi=std(cvnmi)
MFMeasure=mean(cvfm)
SFMeasure=std(cvfm)
save([pwd '/' FolderName{1} '/Regular/Results/all_results.mat'])
