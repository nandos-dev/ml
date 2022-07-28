%multiple genes can targe the same gene 
%exploring the excel. 
%keeping the sample numbers are useful to keep the sample 
%labels to detemine which are outliers 
%Data does not seem to be on a log scale data 

%{
Unsure which normalization technique... 
Do not use the portions of the code, and the assignment solution 

%}

%{
Now we want to extract the Target --> indication of prediction 
Convert all of the ALL match to the 0s and all other methods as 1, as an
example
%} 

%{
have the gene labels and target labels
Moving forward going to partition the data. 
    4th fold of the dataset. Extract 18 samples (serve as the test set of
    the samples to be used)
Perform the extraction a total of 4 times, thus a 4-fold validation. 
7-sample extraction, will be used as the 10-fold validation 
1-sample is the leave-one-out extraction 

4-fold more stringent, but the leave-one-out method is expected to perform
better 
%}

%{

cvpartition, process more involved in the partition of the data, order is
artbitary 
    >provide the target matrix 
cv = cvpartition(T,'k',4)
%return variable cv, ans variables and the indicies.
cv.training(1)
%return a logical array of the 1s which are apart of the training set and
the 0 is not contained in the training set. 

%train set may contain both the training and observations sets. 

The reason why use the cvparition is to avoid and imbalances in the data. 
randomly pic the ALL and AML for the test set. 
%}


%{


Xtrain  = X(train,:); %samples of the training is acrossed all the rows 
size(Xtrain)  
    %ans = 54 samples
%recall that this is the first fold of the cross-validation 


%}


%{

Can choose for the SVM, neural networks, decision tress etc.
Use the logistic regression --> will need to change the fitcsvm
Recall ALWAYS NEED to set the training and the test sample. 


How do we account for handeling very large datasets

%}



%{
mdl.predict(Xtest(1,:)) 
mdl.predict(Xtest) 
    %the sample all came back to be the AML, however the model is 
    %not able to predict them correctly. Not perfect, but some results are
    %produced. 


%}


%{

function numerror = svmdemo_golubb99 = 

%returns the error rate of the methods 



%}


%{

Feature Selection, 
%finding correlation between 0s and 1s 

%Adding to the number of features until the selection features 
%no longer improves. 

svmdemo_golub99_trainandtest(Xtrain(:.1),Ttrain, Xtest(:,1),Ttest)
%repeat these for all the genes, and get to number of genes for the 
%best feature. 



%}

%{

progressing and adding to the number of featurs by having the updating 
Finding the criterion value of hte 0.319 


Ifilter(11) 
fprinf('Feature Selection resulted in the following genes:' \n
, strjoin(genenames(Iselected),', ')); 



%}


mdl = fitcsxm(Xtrain, Ttrain, 'KernalFunction', 'rbf');
Ytest = mdl.predict(Xtest); 
%take the form of the mdl object, oop,  x as an object, similar to 
%structures, methods and fucntions. 

%learning the targets is usually slower. 
%{
processing -- Xtrain, Ttrain (target, input of the variables to 
the genes) 

technical  -- data normalization, and feature selection 

%}


%{
corssval make use of the cross valduation 
will  output the numTestSets and TrainSize and TestSize. 

using some features and not all features. 
accuracy is still 65 percent, even after just looking at the top 50 


%}
corrvals = corr(X, strcmp(T, 'ALL'));
[~,I] = sort(abs(corrvals), 'descend');
Ifilter = I(1:50);


Ifilter_selected = sequnetials(@svmdemo_golub99_trainandtest,X (:,Ifilter, T...
    ,'cv',cv, 'options',statset('display','iter'),'direction','forward') 
Iselected = Ifilter( Ifilter_selected ); 

%{
This occurs because of the noise, if we select all the genes, need to have
the geens or features of most relevant importance. 

Why is normalization important? 




%}

X = bsxfun(@minus,X,mean(X)); 
%{
Taking the X-values, need to avoid the X = mean(x) --> do not used 

X = bsxfun -> shift the numbers to the mean of the values, then divide each
sample with the std.dev across the rows of the samples. 

need to run the crossval to calculate the errors and from there can 
calculate the accuracy, ability for the cross-validation 




%}

%Assignment 
bmes_downloadandparsegse_cached('GSE7390')

%{

patient relapse, on the future to group the patients, and make predictions
of the individual will have a relapse. 
    >purpose based on treatment regimen, and prognosis. 
%}

%76 gene signature 

find(strmp(d.RowNames,'219340_s_at')) 
%separate out the names 


%store the names and store-names into the files. 
%breastcancer76genesignature.txt 

strsplit(fileread(''))

%{

create a map, by a map rather than searching each gene name in the
gene-name matrix. 

first method to use is by using the readtable function 

separate out the text using the newline character 

%}

enteries = strsplit( s, newline) 
% now these enteries for the gene-line 
% store into a single entry 
% now each

% entry = enteries{1}
% strsplit(entry, ' ')

%only each the portion 
%search newline 

regexp(s, '\n[\d_a-zA-Z]+ ','match')
%searching for one or more character 
%newline based on digits and characters 


entires = regexp(s, '(?:\n|^)([\d_a-zA-Z]+)','tokens')
%
%return the tokens, parentheses expressions of the code 
%results are now sotred as cells 


%storing the tokens 
[entries{:}]
%storing these as a single cell-array 

%parenthesized 

%searching the artifical test 
%this has some limitation, taking up the RAM, more memory. 
%memory approach, need to be worried in the gb system requirements. 

enteries = regexp([newline s], '\n[\d_a-zA-Z]+ ','tokens')























