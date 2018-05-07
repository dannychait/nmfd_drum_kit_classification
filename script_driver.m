addpath /Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/util

%Sound files to be decomposed
mixDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/mixes/');
mixFileList = getFileNames(mixDirectory ,'wav');
%Directories of isolated snare sound files%
snareDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/close mics/snares/');
snareFileList = getFileNames(snareDirectory ,'wav');
%Directories of isolated kick sound files%
kickDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/close mics/kicks/');
kickFileList = getFileNames(kickDirectory ,'wav');
%Directories of isolated high hat sound files%
hhDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/close mics/hh/');
hhFileList = getFileNames(hhDirectory ,'wav');
%Mat file directory%
closematFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/closematFiles/');
closematFileList = getFileNames(closematFileDirectory ,'mat');
%Mat file directory, mixes%
mixmatFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/mixmatFiles/');
mixmatFileList = getFileNames(mixmatFileDirectory ,'mat');
%TestMat file directory, mixes%
testmixmatFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/testmixmatFiles/');
testmixmatFileList = getFileNames(testmixmatFileDirectory ,'mat');
%TestMat file directory, close mics%
testclosematFileDirectory = ('/Users/dannychait/Documents/MATLAB/Code/NMF-matlab-master/danny_NMFD/audio/testclosematFiles/');
testclosematFileList = getFileNames(testclosematFileDirectory ,'mat');

%% NMFD Decomp and save mat files %%
for i=1:length(mixFileList)
%for i=1:10
    %% parameters %%
    %length of the FFT%
    Nfft = 512;
    %number of time-frequency atoms%
    R = 10;
    % time-length of the atoms%
    T = 120;
    % number of iterations%
    Niter = 10;
    %% preparing the data %%
    disp('reading mix file...')
    %Load the sound file to be decomposed%
    [s,sr] = audioread([mixDirectory mixFileList{i}]);
    %convert to mono by default%
    s = toMono(s);
    %compute the spectrograms of mix file and close mics%
    sp = stft(s,Nfft,hamming(Nfft,'periodic'),Nfft/4);
    V = abs(sp);
    M = size(V,1);
    N = size(V,2);
    %first decomposition%
    [W,H] = NMFD(V,R,T,Niter);
    %post processing%
    initVal.W = W;
    initVal.H = max(H,max(H(:))/10)-max(H(:))/10+0.00001;
    %second decomposition%
    [W,H] = NMFD(V,R,T,10,initVal);
    %separation of the sounds of each component%
    Lambda = cell(R,1);
    for z = 1:R
        Lambda{z} = zeros(M,N);
        for f = 1:M
            v = reshape(W(f,z,:),T,1);
            cv = conv(v,H(z,:));
            Lambda{z}(f,:) = Lambda{z}(f,:) + cv(1:N);
        end
    end

    LambdaTot = zeros(M,N);
    for z = 1:R
        LambdaTot = LambdaTot +Lambda{z};
    end

    for z = 1:R
        xs{z} = istft(sp.*Lambda{z}./LambdaTot,Nfft,hamming(Nfft,'periodic')',Nfft/4);
        %soundsc(xs{z},sr);
    end

    filename = sprintf('%s_%d','nmfdWH',i);

    disp('storing mix file... ')
    save([mixmatFileDirectory filename '.mat'],'W','H','Lambda','xs')
    
    d = length(s)/sr;

    figure 
    imagesc((1:N)/N*d,0:0.05:sr/2000,db(V))
    title('Original Spectrogram')
    xlabel('time (s)')
    ylabel('frequency (kHz)')

    axis xy
    mx = max(caxis);
    caxis([mx-120,mx])
    dyn = caxis;
    
    figure
    subplot(131)
    imagesc((1:T)/N*d,0:0.05:sr/2000,db(reshape(W(:,1,:),M,T)))
    caxis(dyn)
    axis xy
    title('First template')
    xlabel('time (s)')
    ylabel('frequency (kHz)')

    subplot(132)
    imagesc((1:T)/N*d,0:0.05:sr/2000,db(reshape(W(:,2,:),M,T)))
    caxis(dyn)
    axis xy
    title('Second template')
    xlabel('time (s)')
    ylabel('frequency (kHz)')

    subplot(133)
    imagesc((1:T)/N*d,0:0.05:sr/2000,db(reshape(W(:,3,:),M,T)))
    caxis(dyn)
    axis xy
    title('Third template')
    xlabel('time (s)')
    ylabel('frequency (kHz)')

    figure
    subplot(311)
    plot((1:N)/N*d,(H(1,:)))
    title('Activation of the 1st template')
    xlabel('time (s)')

    subplot(312)
    plot((1:N)/N*d,(H(2,:)))
    title('Activation of the 2nd template')
    xlabel('time (s)')
    subplot(313)
    plot((1:N)/N*d,(H(3,:)))
    title('Activation of the 3rd template')
    xlabel('time (s)')

end


%% Saving Close Mic Files
for i=1:length(snareFileList)
%for i=1:10
    
    disp('reading close files...')
    
    [snclose,fs]=audioread([snareDirectory snareFileList{i}]); 
    [kiclose,fs]=audioread([kickDirectory kickFileList{i}]);  
    [hhclose,fs]=audioread([hhDirectory hhFileList{i}]);
    
    disp('extracting from close files...')
    
    snclose = toMono(snclose); % mono
    kiclose = toMono(kiclose); % mono
    hhclose = toMono(hhclose); % mono
    
    filename = sprintf('%s_%d','closemics',i);
    
    disp('storing close files... ')
    save([closematFileDirectory filename '.mat'],'snclose','kiclose','hhclose')
      
end

%% Reload them into Matlab / SVM

% Main Files
[features,labels,HofLabel] = LoadingProcedure(mixmatFileList,mixmatFileDirectory,closematFileList,closematFileDirectory);

% Test Files
[testfeatures,testlabels,testHofLabel] = LoadingProcedure(testmixmatFileList,testmixmatFileDirectory,testclosematFileList,testclosematFileDirectory);

% save H of labels (locations of 'notes' per label)
disp('storing H of label... ')
save('H_of_label.mat','HofLabel')
% save test H of labels (locations of 'notes' per label)
disp('storing test H of label... ')
save('test_H_of_label.mat','testHofLabel')

% save features and labels
disp('storing features and labels... ')
save('features_labels.mat','features','labels')
% save features and labels
disp('storing test features and labels... ')
save('test_features_labels.mat','testfeatures','testlabels')

% load features and labels, put them in model
load features_labels
X = features;
Y = categorical(labels);
classOrder = unique(Y);
rng(1); % For reproducibility

% do this to the 80 left after pulling 20
t = templateSVM('Standardize',1);
CVMdl = fitcecoc(X,Y,'Learners',t,'ClassNames',classOrder);
CMdl = CVMdl;

load test_features_labels
XTest = testfeatures;
YTest = testlabels';

preds = predict(CMdl,XTest);
table(YTest,preds,...
    'VariableNames',{'TrueLabels','PredictedLabels'})

preds = cellstr(preds);

g1 = YTest';		% Known groups
g2 = preds';	% Predicted groups
[C,order] = confusionmat(g1,g2);

%Look at confusion matrix results
C
order