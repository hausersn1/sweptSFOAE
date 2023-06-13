%% Run swept SFOAE using NEL (but outside it)
% changed to include a nel calib file 3/14/23
% change to be log sweep instead of linear 3/14/23

%% Set up
% Get stimulus structure
stim = makeSFOAEstim_new;
addpath('../FPLclick/')

subj = input('Please subject ID:', 's'); % Create subject
stim.subj = subj;

date = datetime('now');
date.Format = 'yyyy_MM_dd';
datetag = string(date);
stim.date = datetag;

% % Get CalibFile from NEL DATA storage
% ExpDataDir = 'C:\NEL\ExpData\';
% cd(ExpDataDir);
% findDir = dir(sprintf('*%s*%s*', datetag, subj));
% calibflag = 1;
% while calibflag == 1
%     calibNum = input('What NEL calib #? ', 's');
%     if isempty(findDir)
%         fprintf(2,'You need to run a NEL calibration first!\n');
%         return
%     elseif length(findDir)~=1
%         fprintf(2,'Multiple Directories. I am confused. \n')
%         return
%     else
%         datadir = [ExpDataDir findDir.name];
%         cd(datadir);
%         calib_file = dir(sprintf('coef_000%s_calib.mat', calibNum));
%         if isempty(calib_file)
%             if calibNum == '999'
%                 noCalibFlag = 1;
%                 calibflag = 0;
%             else
%                 fprintf(2, 'No calib of that number. Try again! \n')
%             end
%         else
%             load(calib_file.name, 'b');
%             stim.b = b;
%             calibflag = 0;
%         end
%     end
% end


earflag = 1;
while earflag == 1
    ear = input('Please enter which year (L or R):', 's');
    switch ear
        case {'L', 'R', 'l', 'r', 'Left', 'Right', 'left', 'right', 'LEFT', 'RIGHT'}
            earname = strcat(ear, 'Ear');
            earflag = 0;
            stim.ear = ear;
        otherwise
            fprintf(2, 'Unrecognized ear type! Try again!');
    end
end



% Get FPL CalibFile
calibflag = 1;
noCalibFlag = 0;
while calibflag == 1
    stim.b.Ph1 = FPL_inv_calib_fir_coeff(subj, ear, '1'); 
    stim.b.Ph2 = FPL_inv_calib_fir_coeff(subj, ear, '2');  
    calibflag = 0; 
end


% Make directory to save results
mainDir = 'C:\Users\Heinz Lab - NEL2\Desktop\OAEs\SFOAE\swept_wCalib\';
paraDir = [mainDir 'Results'];
cd(mainDir);
addpath(genpath(paraDir));
if(~exist(strcat(paraDir,'\',subj),'dir'))
    mkdir(strcat(paraDir,'\',subj));
end
respDir = strcat(paraDir,'\',subj,'\');

% Remind user about ER-10B+
uiwait(warndlg('Set ER-10B+ GAIN to 40 dB','SET ER-10B+ GAIN WARNING','modal'));

%% Start (w/ Delay if needed)
button = input('Do you want a 10 second delay? (Y or N):', 's');
switch button
    case {'Y', 'y', 'yes', 'Yes', 'YES'}
        DELAY_sec=10;
        fprintf(1, '\n%.f seconds until START...\n',DELAY_sec);
        pause(DELAY_sec)
        fprintf(1, '\nWe waited %.f seconds ...\nStarting Stimulation...\n',DELAY_sec);
    otherwise
        fprintf(1, '\nStarting Stimulation...\n');
end



windowdur = 0.5;
SNRcriterion = 6;
maxTrials = 50;
minTrials = 12;
doneWithTrials = 0;
figure;

% Make arrays to store measured mic outputs
ProbeBuffs = zeros(maxTrials, numel(stim.yProbe));
SuppBuffs = zeros(maxTrials, numel(stim.yProbe));
BothBuffs = zeros(maxTrials, numel(stim.yProbe));
flip = -1;

% Initializing TDT
fig_num=99;
GB_ch=1;
FS_tag = 3;
Fs = 48828.125;
[f1RZ,RZ,~]=load_play_circuit_Nel2(FS_tag,fig_num,GB_ch);

%% Loop for presenting stimuli

% variable for live analysis
k = 1;
t = stim.t;
testfreq = [.75, 1, 1.5, 2, 3, 4, 6, 8, 12].* 1000;

if stim.speed < 0
    f1 = stim.fmax;
    f2 = stim.fmin;
else
    f1 = stim.fmin;
    f2 = stim.fmax;
end

if stim.speed < 20
    t_freq = log2(testfreq/f1)/stim.speed + stim.buffdur;
else
    t_freq = (testfreq-f1)/stim.speed + stim.buffdur;
end


while doneWithTrials == 0
    
    % alternate phase of the suppressor
    flip = flip .* -1;
    
    % Do probe only
    dropSupp = 120;
    dropProbe = stim.drop_Probe;
    buffdata = zeros(2, numel(stim.yProbe));
    buffdata(1, :) = stim.yProbe;
    
    % filter data
    buffdata = filter(stim.b.Ph1, 1, buffdata, [], 1);
    
    % Load the 2ch variable data into the RZ6:
    invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', buffdata(1, :));
    invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', buffdata(2, :));
    % Set the delay of the sound
    invoke(RZ, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
    playrecTrigger = 1;
    % Set attenuations
    rc = PAset([0, 0, dropProbe, dropSupp]);
    % Set total length of sample
    RZ6ADdelay = 97; % Samples
    resplength = size(buffdata,2) + RZ6ADdelay; % How many samples to read from OAE buffer
    invoke(RZ, 'SetTagVal', 'nsamps', resplength);
    
    %Start playing from the buffer:
    invoke(RZ, 'SoftTrg', playrecTrigger);
    currindex = invoke(RZ, 'GetTagVal', 'indexin');
    
    while(currindex < resplength)
        currindex=invoke(RZ, 'GetTagVal', 'indexin');
    end
    
    vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
        'F32','F64',1);
    
    % Save data
    if k > stim.ThrowAway
        ProbeBuffs(k - stim.ThrowAway,  :) = vin((RZ6ADdelay + 1):end);
    end
    % Get ready for next trial
    invoke(RZ, 'SoftTrg', 8); % Stop and clear "OAE" buffer
    %Reset the play index to zero:
    invoke(RZ, 'SoftTrg', 5); %Reset Trigger
    
    pause(0.05);
    
    % Do suppressor only
    dropProbe = 120;
    dropSupp = stim.drop_Supp;
    buffdata = zeros(2, numel(stim.ySupp));
    buffdata(2, :) = flip.*stim.ySupp;
    
    % filter data
    buffdata = filter(stim.b.Ph2, 1, buffdata, [], 1);

    % Load the 2ch variable data into the RZ6:
    invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', buffdata(1, :));
    invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', buffdata(2, :));
    % Set the delay of the sound
    invoke(RZ, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
    playrecTrigger = 1;
    % Set attenuations
    rc = PAset([0, 0, dropProbe, dropSupp]);
    % Set total length of sample
    RZ6ADdelay = 97; % Samples
    resplength = size(buffdata,2) + RZ6ADdelay; % How many samples to read from OAE buffer
    invoke(RZ, 'SetTagVal', 'nsamps', resplength);
    
    %Start playing from the buffer:
    invoke(RZ, 'SoftTrg', playrecTrigger);
    currindex = invoke(RZ, 'GetTagVal', 'indexin');
    
    while(currindex < resplength)
        currindex=invoke(RZ, 'GetTagVal', 'indexin');
    end
    
    vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
        'F32','F64',1);
    
    % Save data
    if k > stim.ThrowAway
        SuppBuffs(k - stim.ThrowAway,  :) = vin((RZ6ADdelay + 1):end);
    end
    % Get ready for next trial
    invoke(RZ, 'SoftTrg', 8); % Stop and clear "OAE" buffer
    %Reset the play index to zero:
    invoke(RZ, 'SoftTrg', 5); %Reset Trigger
    
    pause(0.05);
    
    % Do both
    dropProbe = stim.drop_Probe;
    dropSupp = stim.drop_Supp;
    buffdata = zeros(2, numel(stim.yProbe));
    buffdata(1, :) = stim.yProbe;
    buffdata(2, :) = flip.*stim.ySupp;
    
    % filter data
    buffdata(1,:) = filter(stim.b.Ph1, 1, buffdata(1,:), [], 1);
    buffdata(2,:) = filter(stim.b.Ph2, 1, buffdata(2,:), [], 1);
    
    % Load the 2ch variable data into the RZ6:
    invoke(RZ, 'WriteTagVEX', 'datainL', 0, 'F32', buffdata(1, :));
    invoke(RZ, 'WriteTagVEX', 'datainR', 0, 'F32', buffdata(2, :));
    % Set the delay of the sound
    invoke(RZ, 'SetTagVal', 'onsetdel',100); % onset delay is in ms
    playrecTrigger = 1;
    % Set attenuations
    rc = PAset([0, 0, dropProbe, dropSupp]);
    % Set total length of sample
    RZ6ADdelay = 97; % Samples
    resplength = size(buffdata,2) + RZ6ADdelay; % How many samples to read from OAE buffer
    invoke(RZ, 'SetTagVal', 'nsamps', resplength);
    
    %Start playing from the buffer:
    invoke(RZ, 'SoftTrg', playrecTrigger);
    currindex = invoke(RZ, 'GetTagVal', 'indexin');
    
    while(currindex < resplength)
        currindex=invoke(RZ, 'GetTagVal', 'indexin');
    end
    
    vin = invoke(RZ, 'ReadTagVex', 'dataout', 0, resplength,...
        'F32','F64',1);
    
    % Save data
    if k > stim.ThrowAway
        BothBuffs(k - stim.ThrowAway,  :) = vin((RZ6ADdelay + 1):end);
    end
    % Get ready for next trial
    invoke(RZ, 'SoftTrg', 8); % Stop and clear "OAE" buffer
    %Reset the play index to zero:
    invoke(RZ, 'SoftTrg', 5); %Reset Trigger
    
    pause(0.05);
    
    fprintf(1, 'Done with trial %d / %d\n', k,...
        (stim.ThrowAway + stim.Averages));
    
    % test OAE
    if k - stim.ThrowAway >= minTrials
        SFOAEtrials = ProbeBuffs(1:k, :) + SuppBuffs(1:k, :) - BothBuffs(1:k, :);
        SFOAE = median(SFOAEtrials,1);
        coeffs_temp = zeros(length(testfreq), 2);
        coeffs_noise = zeros(length(testfreq), 8);
        for m = 1:length(testfreq)
            win = find( (t > (t_freq(m)-windowdur/2)) & ...
                (t < (t_freq(m)+windowdur/2)));
            taper = hanning(numel(win))';
            
            resp = SFOAE(win) .* taper;
            
            phiProbe_inst = stim.phiProbe_inst;
            model_sf = [cos(2*pi*phiProbe_inst(win)) .* taper;
                -sin(2*pi*phiProbe_inst(win)) .* taper];
            
            model_noise = [cos(2*pi*1.1*phiProbe_inst(win)) .* taper;
                -sin(2*pi*1.1*phiProbe_inst(win)) .* taper;
                cos(2*pi*1.2*phiProbe_inst(win)) .* taper;
                -sin(2*pi*1.2*phiProbe_inst(win)) .* taper;
                cos(2*pi*1.4*phiProbe_inst(win)) .* taper;
                -sin(2*pi*1.4*phiProbe_inst(win)) .* taper;
                cos(2*pi*1.6*phiProbe_inst(win)) .* taper;
                -sin(2*pi*1.6*phiProbe_inst(win)) .* taper];
            
            coeffs_temp(m,:) = model_sf' \ resp';
            coeffs_noise(m,:) = model_noise' \ resp';
        end
        
        % for noise
        noise2 = zeros(length(testfreq),4);
        for i = 1:2:8
            noise2(:,ceil(i/2)) = abs(complex(coeffs_noise(:,i), coeffs_noise(:,i+1)));
        end
        noise = mean(noise2, 2);
        
        oae = abs(complex(coeffs_temp(:,1), coeffs_temp(:,2)));
        
        SNR_temp = db(oae) - db(noise);
        
        hold off;
        plot(testfreq./1000,db(oae.*10000), 'o', 'linew', 2);
        hold on;
        plot(testfreq./1000,db(noise.*10000), 'x', 'linew', 2);
        legend('SFOAE', 'NOISE');
        xlabel('Frequency (Hz)')
        ylabel('Median Amplitude dB')
        set(gca, 'XScale', 'log', 'FontSize', 14)
        xticks([.5, 1, 2, 4, 8, 16])
        xlim([0.5, 16])
        
        
        if SNR_temp(1:8) > SNRcriterion
            doneWithTrials = 1;
        elseif k == maxTrials
            doneWithTrials = 1;
        end
        
    end
    k = k + 1;
    
end
stim.ProbeBuffs = ProbeBuffs(1:k-1,:);
stim.SuppBuffs = SuppBuffs(1:k-1,:);
stim.BothBuffs = BothBuffs(1:k-1,:);

%% Add useful info to structure
mic_sens = 0.05; % mV / Pa
mic_gain = db2mag(40);

P_ref = 20e-6; % * sqrt(2);
DR_onesided = 1;

stim.VoltageToPascal = 1 / (DR_onesided * mic_gain * mic_sens);
stim.PascalToLinearSPL = 1 /  P_ref;


%% Save Measurements
datetag = datestr(clock);
click.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '_';
fname = strcat(respDir,'SFOAE_log_',subj,'_', earname,'_',datetag, '.mat');
save(fname,'stim');
fprintf(1, 'Saved!');
%% Close TDT, ER-10X connections etc. and cleanup
close_play_circuit(f1RZ, RZ);
rmpath('../FPLclick/'); 

