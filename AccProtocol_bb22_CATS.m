%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Verifying accelerometer signals for calls%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1- Raven selection tables must be merged into 1 file so that start time of ALL selections is relative to 
%start time of first file. That's what you would get if you opened all sound files at once and had one unique table
% for all files. 
    % 1.1 If you have one selectin table for each file, get the total number of
    % samples per audio file (relative to duration of the file) and the
    % sampling rate of the audio. Then, add to the start time (s) of all selections in file n,
    % the total duration (total samples/fs) of file n-1. For example if file 1 is 558307352 samples, and it was sampled at 
    %120.000 Hz. You would add to the start time of selections on file 2 558307352/120000. Do that in R, Matlab or Excel.
    
        %you can get the duration of each file on Raven, when you drag the
        %soundfile to the work window or use [BLKS,fs,fn,recdir,id,recn] = d3getwavcues(recdir,prefix,suffix) or  
        %[ct,ref_time,fs,fn,recdir,recn] = d3getcues(recdir,prefix,suffix)  
    
    %1.2 Create a master Selection table by Merging tables into one using R, Matlab or Excel
    %Make new begin and end time colunms subtracting 30 seconds and add 30
    %s respectively from the start time of the calls. The goal is to have
    %60s selections. Don't change any configurations or titles of the file
    %as Raven has to open it again. Delete the original colunms. You can add more or less time if you think is best.  
    
% 2- Open all the deployment's .wav files in Raven (page it as much as necessary)
% and open the master selection table. Save all selections as .wav files -
% check the file name and MAKE SURE YOU HAVE sel<iii> as the start name of
% the files! IF YOU ONLY HAVE <ii>, THE ACC FILES WILL NOT MATCH THE WAV
% FILES!
% (make sure you have space in the drive to save the selections, depending on how many
%you have they will take a lot of space). 
% Name the folder that contains all selections 'Sound'. 
% 
%% 3- Getting accelerometer samples that match your audio selections:
% Make sure the nyquist of the acc matches the frequency range of the calls
% of interest. I suggest always running the full acc fs with no decimation.
        depid = 'bb22_125a' ;
        prefix=depid;
        %df=1;
        recdir='D:\Tag data\CATS\tag_data_raw\MA\2022\bb20220508-04\raw\wavfiles\';
        %import the .mat file for the CATS deployment using Import Data
        infodep = 'bb'; %had to do this manually, make_info wasn't available for me and I couldn't figure out what the next line needed other than 'bb'
        %A=sens_struct(A,fs,infodep,'acc');
        %A.frame = 'tag' ; 
    %% 3.2 Import selection table 
        tabdir='D:\Tag data\FocalAccAnalysis\bb22_125a\Sound\bb22_125a_noise_removed.txt' %table file with complete address.     
        selections=readtable(tabdir, 'Delimiter', 'tab'); %reading table in. Make sure it is a tab delimited text file. 
        start_time=round(table2array((selections(:,4)))); %For this to be right the selection table MUST have Begin Time on the 4th column
        start_time(:,1)=start_time(:,1)-20 %move start time back 20 seconds and end time forward 20 seconds (line 65) to correspond to a 40s wav file clip
        %and EndTime on the 5th. 
        end_time=round(table2array((selections(:,5))));
        end_time(:,1)=end_time(:,1)+20
        timecues=[start_time end_time]; %time cues must be in seconds 
        
        accdir='D:\Tag data\FocalAccAnalysis\bb22_125a\Acc\' %define directory to save files with acc data. Each file will have 3 variables, one for each acc axis. 
 %% 3.3 Crop samples to start_time and end_time.   
    accsr=400
      
   for v=1:length(timecues);
        A_crop=crop_to(Adata, accsr, [timecues(v,1) timecues(v,2)]); %A must be sensor structure! otherwise you will have to add sf
        accx=A_crop(:,1);
        accy=A_crop(:,2);
        accz=A_crop(:,3);
          if v<10
          fnum=strcat(num2str(0), num2str(0), num2str(v));
          end 
            if 10>=v<=100
            fnum=strcat(num2str(0), num2str(v));
            end
                 if v>100
                 fnum=strcat(num2str(v));
                 end
        fprintf(1, 'saving %s\n', ['file',fnum]); 
            if v<10                   
            save([accdir, '\', 'acc', depid, '0', fnum], 'accx', 'accy', 'accz', 'accsr')
            else
             save([accdir, '\', 'acc', depid, fnum], 'accx', 'accy', 'accz', 'accsr')
            end
        clear ('A_crop', 'accx', 'accy', 'accz', 'fnum') ;
   end
   %% 4 - Run batch processor
    % 4.1 Defining directories  
    %Define directory in which your sound samples are - this is the folder you created in step 2.
    myDirAU='D:\Tag data\FocalAccAnalysis\bb22_125a\Sound'
    %Gets all files in myDirAU that are wav files and stores names in myAUFIles
    myAUFiles=dir(fullfile(myDirAU, '*.wav'));
        
    %Define directory in which acc samples are - folder created in step 3.2
    myDirAC=accdir
    %Gets files with .mat extention and stores names and paths into myACFiles
    myACFiles=dir(fullfile(myDirAC, '*.mat')); %check this structure after to make sure acc files are in sequential order (sometimes 100 comes after 10!)   
    
    %Define directory to save images of the spectrograms 
    specrecdir=('D:\Tag data\FocalAccAnalysis\bb22_125a\Spectrogram\');
    
    %% 4.2 Defining parameters 
    %audio spectrogram 
    freqRange=0:500; %Frequency range of the spectrogram 
    windowsize=5046
    overlap=3800
    
    %%SETTING FILTERS: define filter parameters for acc data
    n=5; %order of filter 
    hp=1;%highpass frequency
    ts=0
    %ACC SPECTROGRAM SETTINGS: Defines parameters for all acc spectrograms. 
    ws=128%window size 
    ov=0.75 %overlap fraction

    %% 4.3 Batch processing 
for k=1:length(myAUFiles)
    baseFileNameAU=myAUFiles(k).name;
    fullFileNameAU=fullfile(myDirAU,baseFileNameAU);
    fprintf(1, 'Working on sound file %s\n', fullFileNameAU); 

    %loads audio file. Gets data and sampling frequency
    [C, fsc]=audioread(fullFileNameAU);

    %gets only one channel (most recordings are stereo)
    C=C(:,1);
        
%AUDIO SPECTROGRAM: Defines figure pannel, calculates spectrogram and plots
%it 
    F1=figure 
    ax(1)=subplot(4,1,1);
    [s,f,t]=spectrogram(C,windowsize, overlap, freqRange, fsc, 'yaxis');
    s = 10*log10(abs(s));
    audiospec = imagesc(0+(t/(60*60*24)),f,s);
    title('AUDIO')

%LOADING IN ACC FILES: for each k file in the directory of ACC samples,
%loads .mat file. Note that k is the same for sound and acc. 
    baseFileNameAC=myACFiles(k).name;
    fullFileNameAC=fullfile(myDirAC,baseFileNameAC);
    fprintf(1, 'Working on acc files %s\n', fullFileNameAC); 
    load(fullFileNameAC);

%APPLYING FILTERS : applies butterworth filter with the specified paramters
%to each acc axis
    [b,a] = butter(n,hp/(accsr/2),'high');
                Aa_fx = filter(b,a,accx);
                Aa_fy = filter(b,a,accy);
                Aa_fz = filter(b,a,accz);
                
%NORMALIZE FILTERED DATA: normalizes data for each axis (takes mean over
%std from each value)
    Aa_fnx=((Aa_fx-mean(Aa_fx)/std(Aa_fx))); %normalizing filtered acc from x axis
    Aa_fny=((Aa_fy-mean(Aa_fy)/std(Aa_fy))); %normalizing filtered acc from y axis
    Aa_fnz=((Aa_fz-mean(Aa_fz)/std(Aa_fz))); %normalizing filtered acc from z axis      


    % Accx: calculates and plots spectrogram for x acc axis
    ax(2)=subplot(4,1,2)
    [s,f,t] = spectrogram(Aa_fnx,hamming(ws),floor(ov*ws),ws,accsr);
            s = 10*log10(abs(s)); % Convert to decibels
            
    accelspecx = imagesc(ts+(t/(60*60*24)),f,s);
    title('Accelerometer X')
%colormap(flipud(bone))

    %Accy: calculates and plots spectrogram for y acc axis
    ax(3)=subplot(4,1,3)
    [s,f,t] = spectrogram(Aa_fny,hamming(ws),floor(ov*ws),ws,accsr);
            s = 10*log10(abs(s)); % Convert to decibels
            
    accelspecy = imagesc(ts+(t/(60*60*24)),f,s);
    title('Accelerometer Y')
    %colormap(flipud(bone))

    %Accz calculates and plots spectrogram for z acc axis
    ax(4)=subplot(4,1,4)
    [s,f,t] = spectrogram(Aa_fnz,hamming(ws),floor(ov*ws),ws,accsr);
            s = 10*log10(abs(s)); % Convert to decibels
            
    accelspecz = imagesc(ts+(t/(60*60*24)),f,s);
    title('Accelerometer Z')
    %colormap(flipud(bone)) 
        
    %SAVING
    %defines name to identify figure based on name of the audio file and acc file
    figname=[baseFileNameAU baseFileNameAC '.png'];
    %saves figure
    print([specrecdir, figname], '-dpng', '-r300')
     close(F1)      
end

fprintf(1, 'JOB DONE...Check specrecdir for results! \n')
        
%%    
list=struct2table(myAUFiles);
list=list(:,1)
writetable(list,'lista.txt')  
        
        
        
        
        
        
        