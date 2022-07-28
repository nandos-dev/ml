function numerror = hwmaml_breastcancer_trainandtest(Xtrain,Ttrain,Xtest,Ttest)

%HWMAML_BREASTCANCER_TRAINANDTEST Summary of this function goes here
%   Detailed explanation goes here

model = fitcsvm(Xtrain,Ttrain, 'KernelFunction', 'rbf');
Ytest = model.predict(Xtest);

numcorrect = sum(strcmp(Ytest,Ttest));
numerror = sum(~strcmp(Ytest,Ttest));
accuracy = numcorrect / numel(Ttest);
fprintf('Accuracy of model is %.2f%%\n',(accuracy*100))
fprintf('numerror of model is %.1f\n',(numerror))
end