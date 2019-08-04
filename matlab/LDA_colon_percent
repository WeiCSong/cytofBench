%% Read the Data and Preprocess
clear
FolderName = {'colon'};
percent_train = 80;
percent_test = 100 - percent_train;
mkdir([pwd '/' FolderName{1} '/' num2str(percent_train) '_percent_train'])

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

DataTable.sample = [];
Data = table2array(DataTable);

clear DataTable

% clear NotDebrisSinglets
Data(strcmp('NotDebrisSinglets',Labels),:)=[];
Labels(strcmp('NotDebrisSinglets',Labels))=[];

% Apply arcsinh5 transformation
Data=asinh((Data-1)/5);

%% run LDA Classifier with 5-fold cross-validation

% CVO = cvpartition(Labels,'HoldOut',percent_test/100);
Accuracy = zeros(5,1);
nmi=zeros(5,1);
FMeasure=zeros(5,1);
training_time = zeros(5,1);
testing_time = zeros(5,1);
%CellTypes = {'','','',};
numCellTypes = unique(Labels);
ConfusionMat = zeros(length(numCellTypes));

    
training_size = round(size(Labels) * percent_train / 100);
training_size = training_size(1);
% testing_size = round(size(Labels) * percent_test / 100);

for i = 1:5
    rperm = randperm(130667, training_size);
    Labels_train = logical(zeros(size(Samples)));
    Labels_train(rperm) = 1;
    Labels_test = logical(ones(size(Samples)));
    Labels_test(rperm) = 0;
    
    trIdx = Labels_train;
    teIdx = Labels_test;
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
save([pwd '/' FolderName{1} '/' num2str(percent_train) '_percent_train/ldaResult_colon_' num2str(percent_train) '_percent.mat'],'Result')
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
saveas(gcf, [pwd '/' FolderName{1} '/' num2str(percent_train) '_percent_train/Performance_Evaluation.png']);

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
saveas(gcf, [pwd '/' FolderName{1} '/' num2str(percent_train) '_percent_train/Population_Frequency.png'])

%% Population Frequency scatter plot
%{
X=log(True_Freq*100);
Y=log(Predicted_Freq*100);
figure,scatter(X,Y,50,'filled')
box on, grid on
xlabel('Log(True frequency %)'),ylabel('Log(Predicted frequency %)')
title('colon')
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
saveas(gcf, [pwd '/' FolderName{1} '/' num2str(percent_train) '_percent_train/Population_Frequency_scatter.png']);
disp(['Training time: ' num2str(training_time)])
disp(['Testing time: ' num2str(testing_time)])
disp(['Total time: ' num2str(Total_time)])

Mac=mean(cvAcc)
Sac=std(cvAcc)
Mnmi=mean(cvnmi)
Snmi=std(cvnmi)
MFMeasure=mean(cvfm)
SFMeasure=std(cvfm)
save([pwd '/' FolderName{1} '/' num2str(percent_train) '_percent_train/all_results.mat'])
