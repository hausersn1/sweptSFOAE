try
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    clear 
    
     initializeER10X;
     % initializeER10X_300Hz_Highpass;
    
    % Initializing TDT
    % Specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    % Get stimulus structure; Change to _linear for linear sweep
    stim = Make_SFswept_log;
    
    if abs(stim.speed) < 20
        sweeptype = 'log'; 
    else 
        sweeptype = 'linear';
    end
        
    % Get subject and ear info
    subj = input('Please subject ID:', 's');
    earflag = 1;
    while earflag == 1
        ear = input('Please enter which ear (L or R):', 's');
        switch ear
            case {'L', 'R', 'l', 'r', 'Left', 'Right', 'left', 'right', 'LEFT', 'RIGHT'}
                earname = strcat(ear, 'Ear');
                earflag = 0;
            otherwise
                fprintf(2, 'Unrecognized ear type! Try again!');
        end
    end
    
    % The button section is just so you can start the program, go into the
    % booth and run yourself as the subject
    button = input('Do you want the subject to press a button to proceed? (Y or N):', 's');
    switch button
        case {'Y', 'y', 'yes', 'Yes', 'YES'}
            getResponse(card.RZ);
            fprintf(1, '\nSubject pressed a button...\nStarting Stimulation...\n');
        otherwise
            fprintf(1, '\nStarting Stimulation...\n');
            fprintf(1, '\nApproximate Duration: %d minutes\n', (stim.t(end).*3.*(stim.ThrowAway + stim.Averages)./60)); 
    end
    
    % Make directory to save results if it doesn't already exist
    paraDir = './Results/';
    
    addpath(genpath(paraDir));
    if(~exist(strcat(paraDir,'\',subj),'dir'))
        mkdir(strcat(paraDir,'\',subj));
    end
    respDir = strcat(paraDir,'\',subj,'\');
    
    %% Present SFOAE stimuli one trial at a time
    % Make arrays to store measured mic outputs
    ProbeBuffs = zeros(stim.Averages, numel(stim.yProbe));
    SuppBuffs = zeros(stim.Averages, numel(stim.yProbe));
    BothBuffs = zeros(stim.Averages, numel(stim.yProbe));
    flip = -1; 
    
    for k = 1:(stim.ThrowAway + stim.Averages)
        
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
            ProbeBuffs(k - stim.ThrowAway,  :) = vins;
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
            SuppBuffs(k - stim.ThrowAway,  :) = vins;
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
            BothBuffs(k - stim.ThrowAway,  :) = vins;
        end
        WaitSecs(0.25);
        
        fprintf(1, 'Done with trial %d / %d\n', k,...
            (stim.ThrowAway + stim.Averages));
        
    end
    stim.ProbeBuffs = ProbeBuffs;
    stim.SuppBuffs = SuppBuffs;
    stim.BothBuffs = BothBuffs;
    
    %% Add useful info to structure
    mic_sens = 50e-3; % mV/Pa
    mic_gain = db2mag(gain + 6); % +6 for balanced cable
    P_ref = 20e-6;
    DR_onesided = 1;
    stim.VoltageToPascal = 1 / (DR_onesided * mic_gain * mic_sens);
    stim.PascalToLinearSPL = 1 /  P_ref;
    
    %% Save Measurements
    datetag = datestr(clock);
    click.date = datetag;
    datetag(strfind(datetag,' ')) = '_';
    datetag(strfind(datetag,':')) = '_';
    fname = strcat(respDir,'SFOAE_',sweeptype,'_',subj,earname,'_',datetag, '.mat');
    save(fname,'stim');
    
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

