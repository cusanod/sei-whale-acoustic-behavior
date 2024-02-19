%%Getting call depth
%DAC Feb 2023
%% Load only the file you are working with - naming them all the same keeps the rest of the code much shorter
clear;

% dtagdata = load('bb_115aData.mat');
% dtagdata = load('bb_115eData.mat');
% dtagdata = load('bb_121aData.mat');
 dtagdata = load('bb_121bData.mat');
% dtagdata = load('bb_121dData.mat');
% dtagdata = load('bb_125bData.mat');
% dtagdata = load('bb_125cData.mat'); 
% dtagdata = load('bb_125eData.mat');
% %catsdata = load('20220502_47.mat'); %no calls on this tag
% catsdata = load('bb_125dData.mat');
% catsdata = load('bb_125aData.mat');
% Hist = load('HistogramData.mat');

%load the associated selection table
% filename = 'bb22_115a_times.txt';
% filename = 'bb22_115e_times_adjusted.txt';
% filename = 'bb22_115e_times_surfacings.txt'; %as a check
% filename = 'bb22_121a_times_adjusted.txt';
% filename = 'bb22_121a_times_surfacings.txt'; %as a check
 filename = 'bb22_121b_times_adjusted.txt';
% filename = 'bb22_121b_times_surfacings.txt'; %as a check
% filename = 'bb22_121d_times_adjusted.txt';
% filename = 'bb22_121d_times_surfacings.txt'; %as a check
% filename = 'bb22_125a_times_adjusted.txt';
% filename = 'bb22_125a_times_surfacings.txt'; %as a check 
% filename = 'bb22_125b_times_adjusted.txt';
% filename = 'bb22_125b_times_surfacings.txt'; %as a check
% filename = 'bb22_125c_times.txt'; 
% filename = 'bb22_125c_times_surfacings.txt'; %as a check - do not use adjusted times 
% filename = 'bb22_125d_times_adjusted.txt';
% filename = 'bb22_125d_times_surfacings.txt'; %as a check
% filename = 'bb22_125e_times_adjusted.txt';
% filename = 'bb22_125e_times_surfacings.txt'; %as a check
%filename = load('20220502_47.txt'); %no calls on this tag

opts = delimitedTextImportOptions("NumVariables", 1);
opts.VariableTypes = "double"; opts.DataLines = [1, Inf];
tagtimes = readtable(filename, opts);
tagtimes = table2array(tagtimes);
%% %create dataframe with depth and time in seconds

% resample data to every 1 second 
% for better results, for resampling when decimation factor r is larger than 13, divide r into smaller factors and call decimate twice
% see help on decimate and also http://animaltags.org/doku.php?id=tagwiki:tools:processing:decdc
resampD = decdc(decdc(dtagdata.depth,5),10); %sample freq was 50Hz here
plott(dtagdata.depth,50,0)
resampD = decdc(decdc(catsdata.data.Pressure,5),10); %sample freq was 50Hz here
plott(catsdata.data.Pressure,50,1)

%round times to nearest integer
tagtimes = round(tagtimes); 
%add column of zeros in the beginning of tag times to prep to add depth
zc = zeros(size(tagtimes,1),1);
tagtimes = [zc, tagtimes];
%create matrix with depth and time
pp=(1:length(resampD))'; %creating dummy x variable same length as resampD 
timedepth = [resampD,pp];

%loop to find times in selection table that match times in tag data
for ii = 1:size(tagtimes,1)
   % Find the index location where 2nd columns match
   idxLocation = timedepth(:,2)==tagtimes(ii,2);  
   % If valid, set the depth (1st column) of the new matrix
   if any(idxLocation)
        tagtimes(ii,1) = timedepth(idxLocation,1);
   end
end

%make and export the table as a csv
header = ["depth","time"]
results = array2table(tagtimes,'VariableNames',header);
writetable(results, 'D:\Tag data\DepthRecal\results\CallDepth_125a.csv') %this also needs to be changed every time because MATLAB sucks

%% Plots
%for this particular project, where we are unsure about the reliability of
%the depth measurements, I have been plotting the dive profiles to visually
%inspect whether the surface is at 0 or not
%in the case of 115a for ex, the dive profile doesn't look 'normal' until
%about 9000 seconds into the record. Until then it never reaches 0 (the
%surface), even though there were in fact surfacings (as per the audio)
%possible the tag slipped? or there is some other reason? unsure

%DTAG
dtagdata.sampleFreq = 50; %append frequency to the structure everytime 
auto_orient_res_plt(dtagdata)

%plot just depth profile
time = (1:length(resampD));
time_unit = 's';
plot(time, resampD*1, 'k-', 'lineWidth', 1.2)%make sure the Depth factor is correct
grid on
ylabel('Depth (m)', 'FontSize', 12)
xlabel(['Time (', time_unit, ')'], 'FontSize', 12)

% CATS
depth = table2array(catsdata.data(:,17));
timesec = 1/catsdata.Hzs.datafs : 1/catsdata.Hzs.datafs : length(depth)/catsdata.Hzs.datafs;
timehour = timesec/3600;

ymin = -25; ymax = 2; xmin = 0; xmax = 17;
figure;
%subplot(3,1,1);
plot(timesec , depth*-1);
ylim([ymin ymax]) ; xlim([xmin xmax]);
ylabel('Depth (m)');
xlabel('Time (Hours)');

%% for paper figures
%calls overlaid on plot for 121a
%load .mat file for tag data from first section
resampD = decdc(decdc(dtagdata.depth,5),10); %sample freq was 50Hz here
%import CallDepth_FocalOnly_121a selection table - right click and import
%data
ds = CallDepthFocalOnly121a(CallDepthFocalOnly121a.calltype=='Downsweep',:);
pulse = CallDepthFocalOnly121a(CallDepthFocalOnly121a.calltype=='Pulse',:);
time = (1:length(resampD));
time_unit = 's';

figure
lw = 0.25;
font_size = 12;
plot(time, resampD*1, 'k-', 'lineWidth', lw)%make sure the Depth factor is correct
grid on
ylabel('Depth (m)', 'FontSize', 12)
xlabel(['Time (', time_unit, ')'], 'FontSize', 12)
hold
%color by call type
plot(ds.time,ds.depth,'pentagram',MarkerFaceColor='#4B0055',MarkerEdgeColor='#4B0055',MarkerSize=15)
plot(pulse.time,pulse.depth,'pentagram',MarkerFaceColor='#FDE333',MarkerEdgeColor='#FDE333',MarkerSize=15) %yellow#FDE333

%121b - good for pulses
%load .mat file for tag data from first section
resampD = decdc(decdc(dtagdata.depth,5),10); %sample freq was 50Hz here
%import CallDepth_FocalOnly_121b selection table - right click and import
%data
CallDepthFocalOnly121b.time = CallDepthFocalOnly121b.time + 41499;
CallDepthFocalOnly121b.time = datetime(CallDepthFocalOnly121b.time, 'ConvertFrom', 'epochtime', 'Epoch', 0, 'Format', 'hh:mm');

knock = CallDepthFocalOnly121b(CallDepthFocalOnly121b.calltype=='Knock',:);
dp = CallDepthFocalOnly121b(CallDepthFocalOnly121b.calltype=='Double pulse',:);
pulse = CallDepthFocalOnly121b(CallDepthFocalOnly121b.calltype=='Pulse',:);
time = (1:length(resampD));

resampT = decdc(decdc(dtagdata.timeHour,5),10);
startTimeSeconds = 41499; %seconds from midnight that tag was on whale
timehh = startTimeSeconds:startTimeSeconds + 66527;
ts = datetime(timehh, 'ConvertFrom', 'epochtime', 'Epoch', 0, 'Format', 'hh:mm');
time_unit = 'hh:mm';

figure
lw = 0.25;
font_size = 12;
plot(ts, resampD*1, 'k-', 'lineWidth', lw)%make sure the Depth factor is correct
grid on
ylabel('Depth (m)', 'FontSize', 12)
xlabel(['Time (', time_unit, ')'], 'FontSize', 12)
datetick('x','HH:MM','keeplimits','keepticks')
hold on;
%color by call type
plot(knock.time,knock.depth,'pentagram',MarkerFaceColor='#008A98',MarkerEdgeColor='#008A98',MarkerSize=12) %green
plot(dp.time,dp.depth,'pentagram',MarkerFaceColor='#A6DA42',MarkerEdgeColor='#A6DA42',MarkerSize=12) %light green'#A6DA42'
plot(pulse.time,pulse.depth,'pentagram',MarkerFaceColor='#FDE333',MarkerEdgeColor='#FDE333',MarkerSize=12) %yellow#FDE333
hold on;
x_points = [ts(31522), ts(31522), ts(62902), ts(62902)];  
y_points = [-70, 10, 10, -70];
color = [0, 0, 1];
a = fill(x_points, y_points, color);
a.FaceAlpha = 0.1;
hold on;
x_points_twi = [ts(29602), ts(29602), ts(31522), ts(31522)];  
y_points_twi = [-70, 10, 10, -70];
color_twi = [0, 0, 1];
a_twi = fill(x_points_twi, y_points_twi, color_twi);
a_twi.FaceAlpha = 0.1;
hold on;
x_points_twi2 = [ts(62902), ts(62902), ts(65002), ts(65002)];  
y_points_twi2 = [-70, 10, 10, -70];
color_twi = [0, 0, 1];
a_twi2 = fill(x_points_twi2, y_points_twi2, color_twi);
a_twi2.FaceAlpha = 0.1;



%125e - meh
%load .mat file for tag data from first section
resampD = decdc(decdc(dtagdata.depth,5),10); %sample freq was 50Hz here
%import CallDepth_FocalOnly_125e selection table - right click and import
%data
knock = CallDepthFocalOnly125e(CallDepthFocalOnly125e.calltype=='Knock',:);
ds = CallDepthFocalOnly125e(CallDepthFocalOnly125e.calltype=='Downsweep',:);
dp = CallDepthFocalOnly125e(CallDepthFocalOnly125e.calltype=='Double pulse',:);
pulse = CallDepthFocalOnly125e(CallDepthFocalOnly125e.calltype=='Pulse',:);
time = (1:length(resampD));
time_unit = 's';

figure
lw = 0.5;
font_size = 12;
plot(time, resampD*1, 'k-', 'lineWidth', lw)%make sure the Depth factor is correct
grid on
ylabel('Depth (m)', 'FontSize', 12)
xlabel(['Time (', time_unit, ')'], 'FontSize', 12)
hold
%color by call type
plot(knock.time,knock.depth,'pentagram',MarkerFaceColor='#008A98',MarkerEdgeColor='#008A98',MarkerSize=15) %green
plot(ds.time,ds.depth,'pentagram',MarkerFaceColor='#4B0055',MarkerEdgeColor='#4B0055',MarkerSize=15)
plot(dp.time,dp.depth,'pentagram',MarkerFaceColor='#A6DA42',MarkerEdgeColor='#A6DA42',MarkerSize=15) %light green'#A6DA42'
plot(pulse.time,pulse.depth,'pentagram',MarkerFaceColor='#FDE333',MarkerEdgeColor='#FDE333',MarkerSize=15) %yellow#FDE333

%125c - good for knocks
%load .mat file for tag data from first section
resampD = decdc(decdc(dtagdata.depth,5),10); %sample freq was 50Hz here
%import CallDepth_FocalOnly_125c selection table - right click and import
%data
knock = CallDepthFocalOnly125c(CallDepthFocalOnly125c.calltype=='Knock',:);
dp = CallDepthFocalOnly125c(CallDepthFocalOnly125c.calltype=='Double pulse',:);
pulse = CallDepthFocalOnly125c(CallDepthFocalOnly125c.calltype=='Pulse',:);
time = (1:length(resampD));
time_unit = 's';

figure
lw = 0.5;
font_size = 12;
plot(time, resampD*1, 'k-', 'lineWidth', lw)%make sure the Depth factor is correct
grid on
ylabel('Depth (m)', 'FontSize', 12)
xlabel(['Time (', time_unit, ')'], 'FontSize', 12)
hold
%color by call type
plot(knock.time,knock.depth,'pentagram',MarkerFaceColor='#008A98',MarkerEdgeColor='#008A98',MarkerSize=15) %green
plot(dp.time,dp.depth,'pentagram',MarkerFaceColor='#A6DA42',MarkerEdgeColor='#A6DA42',MarkerSize=15) %light green'#A6DA42'
plot(pulse.time,pulse.depth,'pentagram',MarkerFaceColor='#FDE333',MarkerEdgeColor='#FDE333',MarkerSize=15) %yellow#FDE333

%115e - not good
%load .mat file for tag data from first section
resampD = decdc(decdc(dtagdata.depth,5),10); %sample freq was 50Hz here
%import CallDepth_FocalOnly_115e selection table - right click and import
%data
knock = CallDepthFocalOnly115e(CallDepthFocalOnly115e.calltype=='Knock',:);
ds = CallDepthFocalOnly115e(CallDepthFocalOnly115e.calltype=='Downsweep',:);
dp = CallDepthFocalOnly115e(CallDepthFocalOnly115e.calltype=='Double pulse',:);
pulse = CallDepthFocalOnly115e(CallDepthFocalOnly115e.calltype=='Pulse',:);
time = (1:length(resampD));
time_unit = 's';

figure
lw = 0.5;
font_size = 12;
plot(time, resampD*1, 'k-', 'lineWidth', lw)%make sure the Depth factor is correct
grid on
ylabel('Depth (m)', 'FontSize', 12)
xlabel(['Time (', time_unit, ')'], 'FontSize', 12)
hold
%color by call type
plot(knock.time,knock.depth,'pentagram',MarkerFaceColor='#008A98',MarkerEdgeColor='#008A98',MarkerSize=15) %green
plot(ds.time,ds.depth,'pentagram',MarkerFaceColor='#4B0055',MarkerEdgeColor='#4B0055',MarkerSize=15)
plot(dp.time,dp.depth,'pentagram',MarkerFaceColor='#A6DA42',MarkerEdgeColor='#A6DA42',MarkerSize=15) %light green'#A6DA42'
plot(pulse.time,pulse.depth,'pentagram',MarkerFaceColor='#FDE333',MarkerEdgeColor='#FDE333',MarkerSize=15) %yellow#FDE333

%% to plot all the plots from auto_orient_res_plt.m but with the calls overlaid on the depth profile
% Extract parameters
knock = CallDepthFocalOnly125c(CallDepthFocalOnly125c.calltype=='Knock',:);
dp = CallDepthFocalOnly125c(CallDepthFocalOnly125c.calltype=='Double pulse',:);
pulse = CallDepthFocalOnly125c(CallDepthFocalOnly125c.calltype=='Pulse',:);

sample_freq = dtagdata.sampleFreq;
A_org = dtagdata.accelInterm;
Depth = dtagdata.depth;
roll_niv = dtagdata.roll;
pitch_niv = dtagdata.pitch;
yaw_niv = dtagdata.head;
signal = dtagdata.relPitch;
positive_peaks = dtagdata.upFluke;
negative_peaks = dtagdata.dnFluke;
bool_active = ~dtagdata.inactive;
amplitude_at_peaks = dtagdata.atpeakAmp;
amplitude_filt = dtagdata.filtAmp;
frequency_at_peaks = dtagdata.atpeakFreq;
frequency_filt = dtagdata.filtFreq;

depth_factor = 1; %set to depth_factor only when the graph is original 

% Plot
time = (1:length(A_org))'/(sample_freq);
time_unit = 's';

figure
m_pk = 5;

lw = 1.2;
lw2 = 1.2;

font_size = 12;

i_p = 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
hold on
plot(time, Depth*depth_factor, 'k-', 'lineWidth', lw)%make sure the Depth factor is correct
grid on
ylabel('Depth (m)', 'FontSize', font_size)
hold on
plot(knock.time,knock.depth,'*',Color='#0000FF')
plot(dp.time,dp.depth,'*',Color='#FF0000')
plot(pulse.time,pulse.depth,'*',Color='#FFFF00')

i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
hold on
plot(time, roll_niv, 'k-', 'lineWidth', lw)
plot(time, pitch_niv, 'b-', 'lineWidth', lw)
plot(time, yaw_niv, 'r-', 'lineWidth', lw)
grid on
ylabel('Orientation (deg)', 'FontSize', font_size)
% legend('roll', 'pitch', 'yaw',...
%   'FontSize', font_size, 'Orientation','horizontal')
legend('roll', 'pitch', 'yaw', 'FontSize', font_size)
% legend('boxoff')

i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
plot(time, signal, 'k-', 'lineWidth', lw2)
hold on
plot(time, positive_peaks, '^')
plot(time, negative_peaks, 'v')
plot(time(~bool_active), signal(~bool_active), 'r.')

grid on
ylabel('Pitch_{dp} (deg)', 'FontSize', font_size)
% legend('signal', 'positive peak', 'negative peak', 'inactive',...
%   'FontSize', font_size, 'Orientation','horizontal')
legend('signal', 'positive peak', 'negative peak', 'inactive', 'FontSize', font_size)
% legend('boxoff')

i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
plot(time, frequency_at_peaks, 'bo', 'lineWidth', lw2)
hold on
plot(time, frequency_filt, 'k-', 'lineWidth', lw2)
grid on
ylabel('Fluke Freq. (Hz)', 'FontSize', font_size)
% legend('at-peak', 'filtered', 'FontSize', font_size, ...
%   'Orientation','horizontal')
legend('at-peak', 'filtered', 'FontSize', font_size)
% legend('boxoff')

i_p = i_p + 1;
ax_pk(i_p) = subplot(m_pk,1,i_p);
plot(time, amplitude_at_peaks, 'bx', 'lineWidth', lw2)
hold on
plot(time, amplitude_filt, 'k-', 'lineWidth', lw2)
grid on
xlabel(['Time (', time_unit, ')'], 'FontSize', font_size)
ylabel('Fluke Amp. (deg)', 'FontSize', font_size)
% legend('at-peak', 'filtered', 'FontSize', font_size, ...
%   'Orientation','horizontal')
legend('at-peak', 'filtered', 'FontSize', font_size)
% legend('boxoff')

linkaxes(ax_pk, 'x')
%% Time at depth (for beh paper)

 dtagdata = load('bb_115aData.mat');
% dtagdata = load('bb_115eData.mat');
% dtagdata = load('bb_121aData.mat');
% dtagdata = load('bb_121bData.mat');
% dtagdata = load('bb_121dData.mat');
% dtagdata = load('bb_125bData.mat');
% dtagdata = load('bb_125cData.mat'); 
% dtagdata = load('bb_125eData.mat');
% %catsdata = load('20220502_47.mat'); %no calls on this tag
% catsdata = load('bb_125dData.mat');
% catsdata = load('bb_125aData.mat');

% resample data to every 1 second 
% for better results, for resampling when decimation factor r is larger than 13, divide r into smaller factors and call decimate twice
% see help on decimate and also http://animaltags.org/doku.php?id=tagwiki:tools:processing:decdc
resampD = decdc(decdc(dtagdata.depth,5),10); %sample freq was 50Hz here
resampD = round(resampD); 
% find the duration of the tag deployment
tagdur = length(resampD);
% time into tag record in seconds where the tag was on the whale and where tag fell off
tagon = 0 %changes for every tag
tagoff = 4 %time from start of tag record, not tagon - changes for every tag
% remove beginning of tag record where the tag was not yet on the animal
% plus 5 minutes to ensure they resumed pre-tagging behavior (based on observations of near immediate return to feeding) 
newdepth = resampD(tagon+300:tagoff); %remove first 0.5 and last 0.5s of clip (padded with noise on either side_


%create matrix with depth and time - do we care?
pp=(1:length(resampD))'; %creating dummy x variable same length as resampD 
timedepth = [resampD,pp];


%% Plot histogram

Foraging = Hist.Foraging;
Traveling = Hist.Traveling;
Resting = Hist.Resting;

figure;
hold on; 
h1 = histfit(Foraging,26,'rayleigh'	);
set(h1(1),'facecolor',	"#0072BD", 'FaceAlpha',0.25, "EdgeColor","#FFFFFF"); 
set(h1(2),'color',	"#0072BD");

h2 = histfit(Traveling,76, 'rayleigh');
set(h2(1),'facecolor',		"#77AC30", 'FaceAlpha',0.25, "EdgeColor","#FFFFFF"); 
set(h2(2),'color',"#77AC30");

h3 = histfit(Resting,73,'rayleigh'	);
set(h3(1),'facecolor',	"#D95319" , 'FaceAlpha',0.25, "EdgeColor","#FFFFFF");
set(h3(2),'color',	"#D95319");

labels1 = {'Foraging','','Traveling','','Resting',''};
xlim([0 15]); lgd = legend(labels1);lgd.FontSize = 12;
ylabel({'Number of Data Points';'Sampling Rate 50Hz'}), xlabel('Depth (Meters)')
hold off;

figure;
h4 = pie([length(Foraging) length(Traveling) length(Resting)]); 
labels2 = {'Foraging','Traveling','Resting'};
lgd = legend(labels2);

T = table([mean(Foraging) ; length(Foraging)/180000], ...
    [mean(Traveling);length(Traveling)/180000], [mean(Resting);length(Resting)/180000],...
    'VariableNames',{'Foraging', 'Traveling', 'Resting'}, ...
    'RowNames',{'Average Depth (Meters)','Total time (Hours)'});
disp(T)

