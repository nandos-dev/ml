function [numerrortrain accuracytrain] = hwmaml_breastcancer_trainandtest(Xtrain,Ttrain,Xtest,Ttest)

%HWMAML_BREASTCANCER_TRAINANDTEST Summary of this function goes here
%   Detailed explanation goes here

model = fitcsvm(Xtrain,Ttrain, 'KernelFunction', 'rbf');
Ytest = model.predict(Xtest);

numcorrect = sum(strcmp(Ytest,Ttest));
numerrortrain = sum(~strcmp(Ytest,Ttest));
accuracytrain = numcorrect / numel(Ttest);
end