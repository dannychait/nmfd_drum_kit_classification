function [features,labels,HofLabel] = LoadingProcedure(mixmatFileList,mixmatFileDirectory,closematFileList,closematFileDirectory)

Nfft=512;
R=10;
close_mics=3;
j=1;

for i=1:length(mixmatFileList)   
    
        load([mixmatFileDirectory mixmatFileList{i}])
    
        load([closematFileDirectory closematFileList{i}])    
 
        spsclose = stft(snclose,Nfft,hamming(Nfft,'periodic'),Nfft/4);
        spkclose = stft(kiclose,Nfft,hamming(Nfft,'periodic'),Nfft/4);
        sphclose = stft(hhclose,Nfft,hamming(Nfft,'periodic'),Nfft/4);
 
        closeMicNames = {'kick','snare','high hat'};
        matches = {'kick','snare','high hat'}; %Default Matching 
        
        l = zeros(3,R);
        
        for z = 1:R
            %max(-100, db(X))
            %l(1,z) = sum(sum((max(-100, abs(spkclose)) - max(-100, abs(Lambda{z}))).^2,1),2);
            l(1,z) = sum(sum((max(-100, db(spkclose)) - max(-100, db(Lambda{z}))).^2,1),2);
            %l(2,z) = sum(sum((max(-100, abs(spsclose)) - max(-100, abs(Lambda{z}))).^2,1),2);
            l(2,z) = sum(sum((max(-100, db(spsclose)) - max(-100, db(Lambda{z}))).^2,1),2);
            %l(3,z) = sum(sum((max(-100, abs(sphclose)) - max(-100, abs(Lambda{z}))).^2,1),2);
            l(3,z) = sum(sum((max(-100, db(sphclose)) - max(-100, db(Lambda{z}))).^2,1),2);
        end
        
        for z = 1:R
            [~,ind] = min(l(:,z));
            matches{z} = closeMicNames{ind}; 
        end

        for f = 1:R
            
            features(j,:) = reshape(W(:,:,f), 1, []);
            %save H peaks at same index
            HofLabel{j} = H(f,:);
        
            labels{j} = matches{f};
        
     
            j = j + 1;
        end
      
    end
end
