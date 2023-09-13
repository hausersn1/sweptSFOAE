%% SFOAEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Samantha Hauser, AuD
% Modified from: Hari Bharadwaj, PhD (SNAP Lab)
% Created: November 2021
% Last revision: 6-Sep-2023 (added artifact checking)
%
% References:
%   SNR based endpoint: Abdala, C., Luo, P. & Shera, C.A. Characterizing
%       the Relationship Between Reflection and Distortion Otoacoustic
%       Emissions in Normal-Hearing Adults. JARO 23, 647–664 (2022).
%       https://doi.org/10.1007/s10162-022-00857-z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data Storage

% Change to 0 if not running on SNAPlab computer at Purdue
Purdue_SNAP = 0;

% Measure-General info
info.name = 'SFOAEswept';
info.version = 'Auto_v01';


% Visit info
if Purdue_SNAP == 1
    % Visit info
    if exist('C:\Experiments\Sam\current_visit.mat','file')
        load('C:\Experiments\Sam\current_visit.mat', 'visit')
        ask = questdlg(sprintf('Is this subject %s?', visit.subj.ID), 'Check Subject', 'Yes', 'No', 'No');
    else
        ask = 'No';
    end
    
    if strcmp(ask, 'No')
        cd ..
        startVisit
        cd(info.name)
    end
    
    subj = visit.subj;
    info.room = visit.room;
    info.univ = visit.univ;
    info.researcher = visit.researcher;
    
else
    subj = input('Subject ID:', 's');
    info.room = 'Room Number';
    info.univ = 'Pitt';
    info.researcher = input('Initials of tester','s');
end

% Get ear info
subj.ear = questdlg('Which ear?', 'Ear', 'L', 'R', 'R');

% Get date/time
datetag = datestr(clock);
info.date = datetag;
datetag(strfind(datetag,' ')) = '_';
datetag(strfind(datetag,':')) = '-';

% Make directory to save results
paraDir = 'C:\Experiments\Sam\SFOAEswept\Results\';
respDir = strcat(paraDir,filesep,visit.subj.ID,filesep);
addpath(genpath(paraDir));
if(~exist(respDir,'dir'))
    mkdir(respDir);
end

fname = strcat(respDir, info.name, '_', ...
    subj.ID, '_', subj.ear, '_', datetag, '.mat');

%% Run Test
tic;
try
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    
    initializeER10X;
    %initializeER10X_300Hz_Highpass;
    
    % Initializing TDT and specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    % Get stimulus structure; Change to _linear for linear sweep
    stim = Make_SFswept;
    
    % Set live analysis parameters
    windowdur = stim.windowdur;
    SNRcriterion = stim.SNRcriterion;
    maxTrials = stim.maxTrials;
    minTrials = stim.minTrials;
    
    phiProbe_inst = stim.phiProbe_inst*2*pi;
    t = stim.t;
    testfreq = stim.testfreq;
    npoints = length(testfreq);
    nearfreqs = stim.nearfreqs;
    VtoSPL = stim.VoltageToPascal .* stim.PascalToLinearSPL;
    
    if stim.speed < 0
        f_start = stim.fmax;
        f_end = stim.fmin;
    else
        f_start = stim.fmin;
        f_end = stim.fmax;
    end
    
    if strcmp(stim.scale, 'log')
        t_freq = log2(testfreq/f_start)/stim.speed + stim.buffdur;
    else
        t_freq = (testfreq-f_start)/stim.speed + stim.buffdur;
    end
    
    k = 0;
    doneWithTrials = 0;
    figure;
    
    % Make arrays to store measured mic outputs
    ProbeBuffs = zeros(maxTrials, numel(stim.yProbe));
    SuppBuffs = zeros(maxTrials, numel(stim.yProbe));
    BothBuffs = zeros(maxTrials, numel(stim.yProbe));
    flip = -1;
    
    while doneWithTrials == 0
        k = k+1;
        k_kept = k - stim.ThrowAway;
        
        % alternate phase of the suppressor
        flip = flip .* -1;
        delayComp = 1; % Always
        
        % Do probe only
        dropSupp = 120;
        dropProbe = stim.drop_Probe;
        buffdata = zeros(2, numel(stim.yProbe));
        buffdata(1, :) = stim.yProbe;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            ProbeBuffs(k_kept,  :) = vins;
        end
        
        WaitSecs(0.25);
        
        % Do suppressor only
        dropProbe = 120;
        dropSupp = stim.drop_Supp;
        buffdata = zeros(2, numel(stim.ySupp));
        buffdata(2, :) = flip.*stim.ySupp;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            SuppBuffs(k_kept,  :) = vins;
        end
        
        WaitSecs(0.25);
        
        % Do both
        dropProbe = stim.drop_Probe;
        dropSupp = stim.drop_Supp;
        buffdata = zeros(2, numel(stim.yProbe));
        buffdata(1, :) = stim.yProbe;
        buffdata(2, :) = flip.*stim.ySupp;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            BothBuffs(k_kept,  :) = vins;
        end
        
        WaitSecs(0.25);
        
        if k > stim.ThrowAway
            % Set empty matricies for next steps
            coeffs_resp = zeros(npoints, 2);
            coeffs_noise = zeros(npoints, 8);
            a_temp = zeros(1, npoints);
            b_temp = zeros(1, npoints);
            
            for p = 1:npoints
                win = find( (t > (t_freq(p) - windowdur/2)) & ...
                    (t < (t_freq(p) + windowdur/2)));
                taper = hanning(numel(win))';
                
                probe_resp = ProbeBuffs(k_kept,:) + SuppBuffs(k_kept,:) - BothBuffs(k_kept,:);
                resp_trial = probe_resp(win).* taper;
                
                model_probe = [cos(phiProbe_inst(win)) .* taper;
                    -sin(phiProbe_inst(win)) .* taper];
                model_noise = ...
                    [cos(nearfreqs(1)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(1)*phiProbe_inst(win)) .* taper;
                    cos(nearfreqs(2)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(2)*phiProbe_inst(win)) .* taper;
                    cos(nearfreqs(3)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(3)*phiProbe_inst(win)) .* taper;
                    cos(nearfreqs(4)*phiProbe_inst(win)) .* taper;
                    -sin(nearfreqs(4)*phiProbe_inst(win)) .* taper];
                
                coeffs_resp(p,:) = model_probe' \ resp_trial';
                coeffs_noise(p,:) = model_noise' \ resp_trial';
            end
            
            a_resp(k_kept,:) = coeffs_resp(:, 1);
            b_resp(k_kept,:) = coeffs_resp(:, 2);
            a_n(k_kept,:) = coeffs_noise(:, 1);
            b_n(k_kept,:) = coeffs_noise(:, 2);
            
            % calculate amplitudes
            oae_temp = abs(complex(a_resp, b_resp));
            median_oae = median(oae_temp,1);
            
            noise = zeros(length(testfreq),4);
            for i = 1:2:8
                noise(:,ceil(i/2)) = abs(complex(coeffs_noise(:,i), coeffs_noise(:,i+1)));
            end
            noise = mean(noise, 2);
            
            % median SNR
            SNR_temp = db(median_oae)' - db(noise);
            
            noisy_trials = 0;
            % artifact check
            if k_kept >= stim.minTrials
                std_oae = std(oae_temp,1);
                for r = 1:k_kept
                    for q = 1:npoints
                        if oae_temp(r,q) > median_oae(1,q) + 4*std_oae(1,q)
                            noisy_trials = noisy_trials+1;
                            break;
                        end
                    end
                end
            end
            
            % if SNR is good enough and we've hit the minimum number of
            % trials, then stop.
            if SNR_temp(1:8) >= SNRcriterion
                if k_kept >= minTrials + noisy_trials
                    doneWithTrials = 1;
                end
            elseif k == maxTrials
                doneWithTrials = 1;
            end
            
            pass = (SNR_temp>=SNRcriterion);
            oae_pass = db(median_oae.*VtoSPL);
            oae_fail = db(median_oae.*VtoSPL);
            oae_pass(~pass) = NaN;
            oae_fail(pass) = NaN;
            
            % Plot amplitudes from live analysis
            hold off;
            plot(testfreq./1000,oae_pass, 'o', 'linew', 2, 'color', [0 0.4470 0.7410]);
            hold on;
            plot(testfreq./1000,oae_fail, 'o', 'linew', 2, 'color', 'k'),
            plot(testfreq./1000,db(noise.*VtoSPL), 'x', 'linew', 2, 'color', [0.6350 0.0780 0.1840]);
            hold off;
            legend('SFOAE', '', 'NOISE', 'location', 'southwest');
            title('SFOAE')
            xlabel('Frequency (Hz)')
            ylabel('Median Amplitude (dB SPL)')
            set(gca, 'XScale', 'log', 'FontSize', 14)
            xlim([0.5, 16]);
            xticks([.5, 1, 2, 4, 8, 16]);
            ylim([-30, 30]);
            yticks((-30:6:30))
            grid on;
            drawnow;
            
            fprintf(1, 'Trials run: %d / Noisy Trials: %d \n', k_kept, noisy_trials);
            
        end
    end
    %% Save full data struct
    data.info = info;
    data.stim = stim;
    data.info.subj = subj;
    data.resp.trialsCollected = k_kept;
    data.resp.ProbeBuffs = ProbeBuffs(1:k_kept,:);
    data.resp.SuppBuffs = SuppBuffs(1:k_kept,:);
    data.resp.BothBuffs = BothBuffs(1:k_kept,:);
    data.info.testDur_s = toc;
    
    save(fname,'data');
    
    %% Close TDT, ER-10X connections etc. and cleanup
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    
catch me
    closeER10X;
    closeCard(card);
    rmpath(pcard);
    rethrow(me);
    
end

