function [Inter,pInt]=InterleaveDetection_nii(FilePath,SaveDir)
%% InterleaveDetection
% This M file runs interleave detection on all the files in a given
% directory. and saves the detected interleave for each file name in a
% stats.txt file.  It then returns "Inter", the interleave value detected by this program.
%
% FilePath: the path of the .nii file you wish to load
% (i.e. /home/USER/fMRI_Data/Subject1/REST_BOLD.nii)
% Note that this takes both .nii and .nii.gz files. The extension must be
% included
%
% SaveDir: Just the directory that you want the files to be saved in.
% Entering '/home/USER/fMRI_Data/Subject1/Interleave_Results' will result
% in a directory 


% NOTE: This program will overwrite any of its previous outputs if you use
% it more than once with the same save directory variable.
%
%


%% Clear vars and Set up Directories
clc
close all

%NOTE: These directories are included in the original Zip File.  If you
%move files or folders, these paths will need to be modified.

addpath('OpenNii');


[pathstr,FileName,ext] = fileparts(FilePath);
[path2,fname2,ext2]=fileparts(pathstr);

if ~isempty(find(FileName=='.'))
    
    FileName=FileName(1:find(FileName=='.')-1);
end

try [image]=load_untouch_nii(FilePath);
catch err
    [image]=nii_read_volume(FilePath);
end

V=image.img;
ImStack=double(V);


%% Set up Stats File

 if exist(SaveDir,'dir')==0
     mkdir(SaveDir)
 end

if ~strcmp(SaveDir(end),'/')    % If the save path already ends with a '/'
    SaveDir=[SaveDir,'/'];      % Add a '/' to the end                      
end                             % Otherwise Do not modify, create stats savefile

SaveFile=[SaveDir,FileName,'_Stats.txt']

if exist(SaveFile,'file')==2    % If a stats savefile already exists
    delete (SaveFile);          % Delete it
end

fid=fopen(SaveFile,'w');
fprintf(fid,'File Name \t Interleave \t Significance\n');
myFormat=' %s\t%f\t%E\n';

%% Begin File Analysis
    
    [X Y Z T]=size(ImStack);    % Get dimensions of the image
    
  %% Find Signal Region
    
    % Time-average to create one 3D volume
    mean_movie=mean(ImStack,4);
    
    % Create 1D histogram and search for first peak
    [ftshist, binpos] = hist(reshape(mean_movie,1,[]),100);
    [pks1,locs]=findpeaks(-ftshist);mmin=(min(locs));threshold_min=binpos(mmin);
    
    % Create second 1D histogram and search for Energy above threshold
    [ftshist, binpos] = hist(reshape(mean_movie,1,[]),1000);
    Energy=ftshist.*binpos.^2;sum_Energy=cumsum(Energy)/sum(Energy);
    
    % Remove top .1% energy
    threshold_max=binpos(min(find((sum_Energy>=0.999))));
    
    % Create Mask
    mask=ones(size(mean_movie));
    mask(find(mean_movie<threshold_min))=0;
    mask(find(mean_movie>=threshold_max))=0;
    
    
    
    %% apply the mask to the entire image
    
    for t=1:size(ImStack,4)
        ImStack(:,:,:,t)=ImStack(:,:,:,t).*mask;
    end;
    
    %% Calculate Slice Signal Average for each Z over all T
    
    for nT=1:T      % for each T
        for n=1:Z   % for each Z
            ave(n,nT)=mean(mean(ImStack(:,:,n,nT))); % Average over x and y
        end
    end
    
           
    % Check for any all zero slices which would lead to "NaN" in statistical analysis
    
    m=1;    
    for n=1:Z
        if sum(1.0*(ave(n,:)~=0))==0;   % If the sum of any one slice average over all time = 0
            ls(m)=n;                    % Then it is all zero, store that slice number in ls(m)
            m=m+1;                      % Incrament m only if an all zero slice was detected
        end
    end
    
    if m~=1             % If m was incramented, ther emust be at least one all zero slice
        ave(ls,:)=[];   % Remove all non zero slices completely from the matrix
    end
    
    
    
    %% Calculate Pearson's correlation between each slice signal average
    
    [rcor,pcor]=corr(ave','type','Kendall','rows','pairwise');
    
    
    %% Compare Each slice's correlation with every other slice up to K-1 slices away
    %
    % This is the Temporal Distance Correlation Matrix
    %
    % This creates a matrix where each row represents an absolute distance
    % from any given slice in row 1.  For example, if the first column has
    % slice "1" in the first row, row 2 will be slice "2", row 5 will be
    % slice "5".  If column 2 has slice "3" in the first row, row 2 will be
    % "4", row 5 will be "7".  In other words, in ANY column, the slice in
    % row 5 is ALWAYS 4 slices away from the slice in row 1.  More
    % generally:
    %
    % Row N, where N != 1 is N-1 slices from row 1.
    %
    % Row 1 = M
    % Row 2 = M + 1
    % Row 3 = M + 2
    % ...etc
    %
    
    
    
    %Compare to 10 slices (The last slice can not be examined, K must be
    % 1 greater than the desired slice).
    
    K=11;    
    k=[0:K-1];  % k to index slices m+k
    t=1;
    
    % Compared  slice m with the K slice after
    for n=1:size(rcor,1)-K      % for all slices n that have K slices after them
        for m=1:length(k)       % compare slice n to slice n+k(m)
            rMean(m,t)=rcor(n+k(m),n);
        end
        t=t+1;
    end
    
    % Compared with the K slices before
    k=[0:-1:-K+1];
    for n=K:size(rcor,1)
        for m=1:length(k);
            rMean(m,t)=rcor(n+k(m),n);
        end
        t=t+1;
    end

    %% Take average and difference of all slice correlations
    
    % Average & Std
    rM=mean(rMean,2);       % Average correaltion values for slices 1 to K away from any given slice
    StrMean=std(rMean,0,2); % Standard Deviation of these averages
    StrM=std(rM);           % Standard Deveation of averaged correlation values from 138
    
    % Difference
    drM=diff(rM);           % First Derivative of averaged signal
    StdrM=std(drM);         % Standard Deviation
    
    % Second Difference
    ddrM=zeros(length(drM),1);  % Zero Pad Second Derivative for indexing convenience
    ddrM(2:end)=diff(drM);      % Second Derivative of Averaged signal
    StddrM=std(ddrM);           % Standard Deviation of the Second derivative
    
    %% Statistics and Interleave Detection
    
    % Take derivative of each column in the Temporal Distance Correlation
    % Function
    
    drMean=diff(rMean,1);           % Difference Accross Columns 
    ddrMean=zeros(size(drMean));    % Zero Pad Second Differnce For Indexing Convenience
    ddrMean(2:end,:)=diff(drMean);  % Second Difference Across Columns
    StddrMean=std(ddrMean,0,2);     % Standard Deviations of each row of the second difference
    
    StdrMean=std(drMean,0,2);       % Standard Deviations of each row of the first difference

    
    % Find local Maxima and Minima for the original averaged temporal
    % distance correlation matrix and its second derivative, respectively.
    
    [pks1,idx1]=findpeaks(-1.*ddrM);    % Find local min of the second derivative (Invert function to do so)
    [pks2,idx2]=findpeaks(rM(2:end));   % Find local mad of the original function (Exclude the first value as it is always "1")
    pks1=pks1*-1;                       % Restore original values to local minima by un-inverting

    
    % Examine The second derivative
    % ---------------------------------------------------------------------
    
for n=1:length(idx1)    % For each local minimum found
    
    %
    % IF
    %
    % CONDITION 1&2:
    % the index of the peak is at least one less than the length of the
    % array, AND the index of the peak is at least one greater than 1
    % (If the peak is not the first or last value of the array)
    %
    % AND
    %
    % CONDITION 2:
    % the magnitude of the peak is further away from the mean than a
    % given threshold (1.25x larger)
    %
    %THEN
    %
    % take the geometric average of the significance (P-value) of this peak
    % compared to its left and right neighbors using ANOVA
    %
        
    if (idx1(n)+1<=length(ddrM))&&(idx1(n)-1>=1)&&(abs(pks1(n)-mean(ddrM(3:end)))>abs(mean(ddrM(3:end)).*1.25))
        ddrMin(n)=geomean([anova1(ddrMean(idx1(n)-1:idx1(n),:)',[],'off'),anova1(ddrMean(idx1(n):idx1(n)+1,:)',[],'off')]);
   
    %
    % ELSE IF
    %
    % CONDITION 1:
    % the peak is the last value in the array
    %
    % AND
    %
    % CONDITION 2:
    % the magnitude of the peak is further away from the mean than a
    % given threshold (1.25x larger)
    %
    % THEN
    %
    % Only compare it to the Left neighbor (The only one it has)
    %
    
    elseif (idx1(n)+1>length(ddrM))&&(abs(pks1(n)-mean(ddrM(3:end)))>mean(ddrM(3:end)).*1.25)
        ddrMin(n)=anova1(ddrMean(idx1(n)-1:idx1(n),:)',[],'off'); 
        
    %
    % ELSE
    % It must be the first value.  Due to the fact that every value in the
    % first row of the original matrix all being one, these statistics are
    % useless.  Ignore this peak.
    %
    
    else
        ddrMin(n)=1;
    end
end

[pDif,Int]=min(ddrMin);     % Take the minimum (most significant) value for the interleave value (as determined by the second derivative)
IntD=idx1(Int)-1;           % Correct for shifts in indexing due to derivatives and padding.




 
    % Examine The Original Averaged Temporal Correlation Matrix
    % ---------------------------------------------------------------------

idx2(end+1)=1;  % The interleave may be '1', which wouldn't be picked up as a local peak, so include that value


for n=1:length(idx2)    % for each local maximum found
        
    %
    % IF
    %
    % CONDITION 1&2:
    % the index of the peak is at least one less than the length of the
    % array, AND the index of the peak is at least one greater than 1
    % (If the peak is not the first or last value of the array)
    %
    %THEN
    %
    % take the geometric average of the significance (P-value) of this peak
    % compared to its left and right neighbors using ANOVA
    %
    
    if (idx2(n)+2<=length(rM))&&(idx2(n)-1>=1)        
        rMin(n)=geomean([anova1(rMean(idx2(n):idx2(n)+1,:)',[],'off'),anova1(rMean(idx2(n)+1:idx2(n)+2,:)',[],'off')]);
        
    %
    % ELSE IF
    %
    % CONDITION 1:
    % the peak is the last value in the array
    %
    % THEN
    %
    % Only compare it to the Left neighbor (The only one it has)
    %
        
    elseif idx2(n)+2>length(rM)
        rMin(n)=anova1(rMean(idx2(n):idx2(n)+1,:)',[],'off');
        
    %
    % ELSE
    % It must be the first value.  Due to the fact that every value in the
    % first row of the original matrix all being one, the left most
    % neighbor will unfairly skew the results.  Only compare to the right
    % neighbor.
    %   
        
    else
        rMin(n)=anova1(rMean(idx2(n)+1:idx2(n)+2,:)',[],'off');
    end
end



[p,Int]=min(rMin);      % Take the minimum (most significant) value for the interleave value (as determined by the Original Matrix)
Interleave=idx2(Int);

    %% Compare Interleave Values and Decide on true interleave Value
    %----------------------------------------------------------------------


    % If the same interleave value is predicted by both the 2nd derivative
    % and the original funcion, use that value.
    if Interleave==IntD
        Inter=IntD;
        pInt=min([p,pDif]);
        
        % Check to see if pDif (from 2nd derivative) is weak (relative to a
        % common significance value from this type of data.  Typical values
        % are < 1e-8
        
    elseif (log10(pDif)>=-3)
        
        % If the original function's predicted interleave value is
        % sufficiently more significant (1.25x), then use the value from
        % the original temporal correlation matrix.
        
        if (p*1.25<pDif)
            Inter=Interleave;
            pInt=p;
            
        % Otherwise, stick with the interleave value from the second
        % derivative of the Temporal Correlation Matrix     
        
        else
            Inter=IntD;
            pInt=pDif;
        end
        
        % Even if the significance of the second derivative value is high
        % (from 310, p < 1e-3), compare it to the the significance of the
        % value from the original temporal correlation matrix.  If the
        % significance from the orignial function is more than 1.5x the
        % order of the derivative value, use the interleave from the
        % original function
        
    elseif (abs(log10(pDif).*1.5)+log10(p)<=0)
        Inter=Interleave;
        pInt=p;
    else
        Inter=IntD;
        pInt=pDif;
    end


%% Plot and Save Figures for Analysis


h=figure;

subplot(2,6,1:3); imagesc(rcor); title(FileName);
subplot(2,6,4:6); imagesc(rMean); title('rMean');
subplot(2,6,7:8); errorbar(rM,StrMean); title('mean of rMean'); %xlabel(['Int=',num2str(Interleave),' p=',num2str(p)]);
subplot(2,6,9:10); errorbar(drM,StdrMean); title('mean of rMean 1st Derivative'); xlabel(['INT= ',num2str(Inter), ' p= ',num2str(pInt)]);
subplot(2,6,11:12); errorbar(ddrM,StddrMean); title('mean of rMean 2nd derivative');% xlabel(['Int=',num2str(IntD),' p=',num2str(pDif)]);

print(h,'-dpng',[SaveDir,FileName,'.png']);
close all

    
    
    %% Setup File Write
    
    fid=fopen(SaveFile,'a');
    fprintf(fid,myFormat,FileName,Inter,pInt);
    fclose(fid);
    
    
end
