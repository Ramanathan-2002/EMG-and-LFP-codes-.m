# EMG-and-LFP-codes-.m
%% In the project focus on the animal models EMG-(Electromyography)  and  LFP- (Local Field Potential ) . It is a matlab based Neural signal analysis 
%% EMG_CHN_32 and LFP_CHN_06

clc;
clear all;


%import the EMG_CHANNEL_DATA 
[Timestamps_EMG, ChannelNumbers_EMG, SampleFrequencies_EMG, NumberOfvaildSamples_EMG, Samples_EMG,Header_EMG ] =  Nlx2MatCSC('path of the file',[1 1 1 1 1],1 ,1 ,[]);

%downsampling the EMG channel data for 1khz
Samples_inVolts = Samples_EMG *   0.000000045776367187500001;
samples_inmVolts = Samples_inVolts * 1000;
EMG_Chn32 = reshape(samples_inmVolts, 1, []);
EMG_Chn32_DS = downsample(EMG_Chn32, 32);

%window length 
wl = 250;
[yupper_EMG, ~] = envelope(EMG_Chn32_DS, wl, 'rms');

% Using Threshold value to Cut the sleep data only 
SampleFrequencies_EMG = 32000/32; 
epoch_length = 3 * SampleFrequencies_EMG; % 3-second epoch
threshold_EMG = 0.15; 
consecutive_wake_threshold = 2; 

%store the sleep only 3sec epoch EMG data  
num_epochs = floor(length(yupper_EMG) / epoch_length); 
epochs_states = strings(num_epochs, 1);  
sleep_epoch_data_only_EMG = []; 
consecutive_wake_count = 0; 

%Elapsed time For tic and toc
tic

            %loop for sleep and wake detection
            for i = 1:num_epochs
                start_idx = (i - 1) * epoch_length + 1;
                end_idx = min(start_idx + epoch_length - 1, length(yupper_EMG));
                epoch_data_RMS = yupper_EMG(start_idx:end_idx);       
                

               %EMG in RMS value of the Emg 
                epoch_rms = rms(epoch_data_RMS);
               
                %if 2 consecutive wake count delete the wake state
                if epoch_rms > threshold_EMG
                    state = 'Awake';
                    consecutive_wake_count = consecutive_wake_count + 1;
                else
                    state = 'Sleep';
                    consecutive_wake_count = 0; 
                    sleep_epoch_data_only_EMG = [sleep_epoch_data_only_EMG , epoch_data_RMS]; 
                end
                 
                % Take sleep only and Print statement for check
                epochs_states(i) = state;
                
                if consecutive_wake_count > consecutive_wake_threshold
                    fprintf('More than %d consecutive Wake epochs detected. Skipping Wake epochs.\n', consecutive_wake_threshold);
                    continue; 
                end
           
               % fprintf('Epoch %d : %s (RMS = %.4f)\n', i, state, epoch_rms);

            end
toc




                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LFP_DATA_PROCESS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                    

                                                                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                                                                                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load LFP_CHANNEL_06 

%import the LFP_CHANNELS_DATA

[Timestamps_LFP, ChannelNumbers_LFP, SampleFrequencies_LFP , NumberOfvaildsamples_LFP , Samples_LFP ,Header_LFP] = Nlx2MatCSC('path of the file',[1 1 1 1 1] , 1 , 1 ,[]);


    % 
    %     %% Load LFP_CHANNEL_07
    % 
     % [Timestamps_LFP , ChannelNumbers_LFP , SampleFrequencies_LFP , NumberOfvaildsamples_LFP , Samples_LFP , Header_LFP ] = Nlx2MatCSC('path of the file', [1 1 1 1 1] , 1 , 1 ,[]);
    % 
    %     %% Load LFP_CHANNEL_22
    % 
     %[Timestamps_LFP , ChannelNumbers_LFP , SampleFrequencies_LFP , NumberOfvaildsamples_LFP , Samples_LFP , Header_LFP ] = Nlx2MatCSC('path of the file', [1 1 1 1 1] ,1 ,1 , []);
    % 
    %     %% Load LFP_CHANNEL_23
    % 
     % [Timestamps_LFP , ChannelNumbers_LFP , SampleFrequencies_LFP , NumberOfvaildsamples_LFP ,  Samples_LFP , Header_LFP] = Nlx2MatCSC('path of the file', [1 1 1 1 1 ] , 1 , 1 , []);



% Downsample for the LFP_Data and  Set window length = 250
Samples_inVolts = Samples_LFP *  0.000000045776367187500001;
Samples_inmVolts = Samples_inVolts * 1000;
LFP_Chn07 = reshape (Samples_inmVolts , 1 , []);
LFP_Chn07_DS = downsample(LFP_Chn07 , 32);
wl = 250;
[yupper_LFP,~] = envelope(LFP_Chn07_DS,wl,'rms');

%SampleFrquencies of LFP_Data
SampleFrequencies_LFP = 32000/32;
threshold_LFP = 0.3;

%store the sleep data
sleep_epoch_data_only_LFP = [];
sleep_epoch_data_only_rms_value_LFP = [];
sleep_epoch_RMS_LFP_data_use_threshold= [];

epoch_length = 3* SampleFrequencies_LFP;
num_epochs = floor(length(yupper_LFP) / epoch_length);
epochs_states_LFP = strings(num_epochs, 1);

consecutive_check = 0;
consecutive_wake_threshold_LFP_RMS = 1;


            
  tic   
  
           %loop for clssify the  number of epochs  yupper LFP data
            for i = 1:num_epochs
                start_idx = (i - 1) * epoch_length + 1 ;
                end_idx = min(start_idx  + epoch_length - 1 ,length(yupper_LFP));
                epoch_data_LFP = yupper_LFP(start_idx : end_idx);   

                % EMG-Sleep use to same in LPF Sleep identify
                 if strcmp(epochs_states(i), 'Sleep')
                    sleep_epoch_data_only_LFP = [sleep_epoch_data_only_LFP ,epoch_data_LFP];

                    rms_value_LFP_sleep = sqrt(mean(epoch_data_LFP.^2));
                    sleep_epoch_data_only_rms_value_LFP = [sleep_epoch_data_only_rms_value_LFP , rms_value_LFP_sleep];

                  % Use 3sec epoch in LFP_RMS_Threshold 
                    if rms_value_LFP_sleep > threshold_LFP
                        State = 'Awake';
                        consecutive_check = consecutive_check + 1;
                    else
                        State = 'Sleep';
                        consecutive_check = 0;
                        sleep_epoch_RMS_LFP_data_use_threshold = [sleep_epoch_RMS_LFP_data_use_threshold, epoch_data_LFP];
                    end

                    epochs_states_LFP(i) = State;                  
                    if consecutive_check > consecutive_wake_threshold_LFP_RMS
                       fprintf("single %d consecutive wake epoch detected . skipped Wake epochs.\n " , consecutive_wake_threshold_LFP_RMS);
                        continue;
                    end
                       % fprintf('epoch %d :%s (Thresold = %4f)\n', i,State , rms_value);

                end

            end 

toc         


%time axis for the EMG DATA and LFP Data  plotting
time_axis_EMG = (1:length(yupper_EMG)) / SampleFrequencies_EMG; 
filtered_time_EMG = (1:length(sleep_epoch_data_only_EMG)) / SampleFrequencies_EMG; 

time_axis_LFP = (1:length(yupper_LFP)) / SampleFrequencies_LFP;
Filtered_time_LFP = (1:length(sleep_epoch_data_only_LFP)) / SampleFrequencies_LFP;

time_axis_LFP_rms = (1:length(yupper_LFP)) / SampleFrequencies_LFP ;
Filtered_time_LFP_rms = (1:length(sleep_epoch_data_only_rms_value_LFP)) *epoch_length / SampleFrequencies_LFP;

time_axis_LFP_RMS = (1:length(yupper_LFP)/ SampleFrequencies_LFP);
Filtered_RMS_LFP = (1:length(sleep_epoch_RMS_LFP_data_use_threshold))*epoch_length / SampleFrequencies_LFP;



%% Parameter Compute DELTA AND THETA

% Sampling rate and epoch setup
SampleFrequencies_LFP = 32000 / 32;  % 1kHz after downsampling
epoch_length = 3; % seconds
samples_per_epoch = epoch_length * SampleFrequencies_LFP;
LFP_Signal = sleep_epoch_RMS_LFP_data_use_threshold; % Use your LFP data


% Normalize LFP signal
Normalized_LFP_data = (LFP_Signal / max(LFP_Signal));
Num_epochs = floor(length(Normalized_LFP_data) / samples_per_epoch);

% Define frequency bands
Delta_band = [0.5 4];
Theta_band = [4 8]; 
Sigma_band = [9 17];

% Design bandpass filters
d_delta = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', Delta_band(1), 'HalfPowerFrequency2', Delta_band(2), ...
    'SampleRate', SampleFrequencies_LFP, 'DesignMethod', 'butter');

d_theta = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', Theta_band(1), 'HalfPowerFrequency2', Theta_band(2), ...
    'SampleRate', SampleFrequencies_LFP, 'DesignMethod', 'butter');

d_sigma = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', Sigma_band(1), 'HalfPowerFrequency2', Sigma_band(2), ...
    'SampleRate', SampleFrequencies_LFP, 'DesignMethod', 'butter');

% Initialize power arrays
Delta_Power = zeros(1, Num_epochs);
Theta_Power = zeros(1, Num_epochs);
Sigma_Power = zeros(1, Num_epochs);


for i = 1:Num_epochs
    epoch_data_temp = LFP_Signal((i-1)*samples_per_epoch + 1 : i*samples_per_epoch);
    epoch_data_RMS = epoch_data_temp ./ max(epoch_data_temp); 
    LFP_Delta = filtfilt(d_delta, epoch_data_RMS);
    LFP_Theta = filtfilt(d_theta, epoch_data_RMS);
    LFP_Sigma = filtfilt(d_sigma, epoch_data_RMS);

    Delta_Power(i) = mean(LFP_Delta .^ 2);
    Theta_Power(i) = mean(LFP_Theta .^ 2);
    Sigma_Power(i) = mean(LFP_Sigma .^ 2);
end


DTR = Delta_Power ./ Theta_Power ; 
TDR = Theta_Power ./ Delta_Power ;
SDR = Sigma_Power ./ Delta_Power ;


DTR = (DTR - min(DTR)) / (max(DTR) - min(DTR));  
threshold = mean(DTR) + 1 * std(DTR);  

sleep_stage = strings(1,Num_epochs);
for i = 1:Num_epochs 
    if DTR(i) > 1
        sleep_stage(i) = "REM";
    else
        sleep_stage(i) = "NREM ";
    end 

end 

sleep_stage = strings(1, Num_epochs);
sleep_stage(DTR > threshold) = 'REM';
sleep_stage(DTR <= threshold) = 'NREM';


sleep_numeric = double(DTR > threshold);
time_axis = (1:Num_epochs)*epoch_length;


%Calculating Time Duration for REM Sleep and NREM Sleep 
num_of_epoch_rem_sleep = sum(sleep_numeric == 1);
num_of_epoch_nrem_sleep = sum(sleep_numeric == 0);

rem_sleep_duration = num_of_epoch_rem_sleep * epoch_length/60 ;
nrem_sleep_duration = num_of_epoch_nrem_sleep * epoch_length/60;

fprintf('Total duration of rem sleep: %.2f minute\n', rem_sleep_duration);
fprintf('Total duration of nrem sleep: %.2f minute\n',nrem_sleep_duration);

% Time Axis for Delta and Theta 
time_axis_of_delta_power = (1:Num_epochs ) * epoch_length;
time_axis_of_theta_power = (1:Num_epochs) * epoch_length;
time_axis_of_Delta_Theta_ratio = (1:Num_epochs) * epoch_length;
time_axis_of_sigma_power = (1:Num_epochs) * epoch_length;

%% DELTA AND THETA PLOTS 

%Delta_power and Theta power 
subplot(3,1,1);
plot(time_axis_of_delta_power , Delta_Power, 'b', 'LineWidth', 1.5);
hold on;
plot(time_axis_of_theta_power, Theta_Power, 'r', 'LineWidth', 1.5);
legend('Delta Power', 'Theta Power');
title('Delta & Theta Power ');
xlabel('Epochs');
ylabel('Power');
hold off;

%Delta_Theta_Ratio
subplot(3,1,2);
plot(time_axis_of_Delta_Theta_ratio , DTR );
title('Delta and Theta Ratio');
xlabel('Epochs');
ylabel('Ratio');


% Sigma Power
subplot(3,1,3);
bar(time_axis_of_sigma_power, Sigma_Power,'b');
title('Sigma Power (9-17 Hz)');
xlabel('Epochs');
ylabel('Power');


%plots For label REM and NREM 
figure;
hold on;
plot(time_axis, sleep_numeric,'b','LineWidth',1.5); 

for i = 1:length(sleep_numeric)
    if sleep_numeric(i) == 1
        stairs(time_axis(i), 1, 'r.', 'MarkerSize', 10); 
    end
end

yticks([0 1]);
yticklabels({'NREM', 'REM'});
xlabel('Time (sec)');
ylabel('Sleep Stage');
title('Sleep Stages: REM vs NREM');
xlim([0 max(time_axis)]);
 ylim([0 1.5]);
hold off;



%% Emg plots

figure;
 subplot(4,1,1);
 plot(EMG_Chn32);
 xlabel('Time Duration');
 ylabel('Samples-inmVolts');
 title('EMG RAW DATA CHN-32');

 subplot(4,1,2);
 plot(EMG_Chn32_DS);
 xlabel('Time Duration');
 ylabel("Samples-inmVolts");
 title('EMG DOWNSAMPLE DATA');

 subplot(4,1,3);
 plot(yupper_EMG);
 xlabel('Time Duration');
 ylabel('Samples-inmVolts');
 title('EMG ENVELOPE DATA');

 subplot(4,1,4);
 plot(sleep_epoch_data_only_EMG);
 xlabel('Time Duration');
 ylabel('Samples-inmVolts');
 title('EMG SLEEP DATA ONLY');


 %% LFP plots 

  figure;
   subplot(4,1,1);
   plot(LFP_Chn07);
   xlabel('Time Duration');
   ylabel('Samples-inmVolts');
   title('LFP RAW DATA CHN-06 ');

   subplot(4,1,2);
   plot(LFP_Chn07_DS);
   xlabel('Time Duration');
   ylabel('Samples-inmVolts');
   title('LFP DOWNSAMPLE DATA');

   subplot(4,1,3);
   plot(yupper_LFP);
   xlabel('Time Duration');
   ylabel('Samples-inmVolts');
   title( 'LFP ENVELOPED DATA');

   subplot(4,1,4);
   plot(sleep_epoch_RMS_LFP_data_use_threshold);
   xlabel('Time Duration');
   ylabel('Samples-inmVolts');
   title('SLEEP DATA ONLY FROM LFP ');







  
