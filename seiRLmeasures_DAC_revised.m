% Program modified from Peter Madsen's batchmag.m script. Written to use Mark Johnson's script for postemph.m to correct or compensate for the highpass filter on the Dtag in measurements of the signals. 
%All signals have been decimated to 24kHz sampling rate to eliminate high frequency noise sources.

%Compiled to make RL measurements from sei whale tags by DAC 2022.

%Define directory in which your sound clips are
myDirAU='D:\SeiWhales\Tag data\ReceivedLevel\bb22_125e\RLclipsNopad\knocks'

%Gets all files in myDirAU that are wav files and stores names in myAUFIles
s=dir(fullfile(myDirAU, '*.wav'));

% Enter the correct calibration value here for individual Dtags plus any gain. Otherwise use the 178 value from Holt et al. 2017 calibration of D3's.
cal= -178 %- 12.1;
%ONLY ADD GAIN FOR 115a and 115e - the difference between those deployments and the rest is 12.1! 
%create empty data frames
magsrms=zeros(length(s),1);
ptp=zeros(length(s),1);

for k=1:numel(s)
    baseFileNameAU=s(k).name;
    fullFileNameAU=fullfile(myDirAU,baseFileNameAU);
    fprintf(1, 'Working on sound file %s\n', fullFileNameAU); 

    % loads audio file. Gets data and sampling frequency
    [x, fs]=audioread(fullFileNameAU);

    %postemph undoes the effects of the highpass filter from the Dtags
    %postemph does not work with almost all knocks! they are too short! don't use for those
    y=x; %if not using postemph
    %[y,fs] = d3postemph(x,fs); %do not use with CATS
    %y = decimate(y,3); %decimation - NEED FOR CATS ONLY to make all clips 24kHz sampling rate

    % Calculate energy
    CE=cumsum(y.^2); % accumulate energy in window
    ce= sum(y.^2); %sum energy
    ce=CE/ce; %normalized energy
    k1= min(find(ce>0.05*(max(ce)))) ; %getting value that is the 5th percentile
    k2= min(find(ce>0.95*(max(ce)))) ; %getting value that is the 95th percentile
    y=y(k1:k2); %window containing 90% energy

    %RMS 
    % Calculate the RMS value of the acoustic pressure
    SPL = 20*log10(sqrt(mean(y.^2))) - cal;
    magsrms(k)=SPL;

    %peak to peak SPL (for impulsive sounds, knocks and pulses)
    pktopk = 20*log10(max(y)-min(y)) - cal;
    ptp(k)=pktopk;
end

%note that these matrices are saved with clips that might be out of order! Need
%to take care to line up the sound files used with the proper values - see
%command window to see what clips were processed in what order
writematrix(magsrms,'D:\SeiWhales\Tag data\ReceivedLevel\results\bb22_125eRL_rms_knocks.txt')
writematrix(ptp,'D:\SeiWhales\Tag data\ReceivedLevel\results\bb22_125eRL_ptp_knocks.txt')


clear
