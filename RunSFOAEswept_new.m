try
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    initializeER10X;
    
    % Initializing TDT
    % Specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    
    % Get stimulus structure
    stim = makeSFOAEstim_new;
    
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
    ProbeBuffs = zeros(stim.Averages, numel(stim.t));
    SuppBuffs = zeros(stim.Averages, numel(stim.t));
    BothBuffs = zeros(stim.Averages, numel(stim.t));
    
    for k = 1:(stim.ThrowAway + stim.Averages)
        delayComp = 1; % Always
        % Do probe only
        dropSupp = 120;
        dropProbe = stim.drop_Probe;
        buffdata = zeros(2, numel(stim.t));
        buffdata(1, :) = stim.yProbe;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        
        if k > stim.ThrowAway
            ProbeBuffs(k - stim.ThrowAway,  :) = vins;
        end
        
        % Do suppressor only
        dropProbe = 120;
        dropSupp = stim.drop_Supp;
        buffdata = zeros(2, numel(stim.t));
        buffdata(2, :) = stim.ySupp;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            SuppBuffs(k - stim.ThrowAway,  :) = vins;
        end
        
        % Do both
        dropProbe = stim.drop_Probe;
        dropSupp = stim.drop_Supp;
        buffdata = zeros(2, numel(stim.t));
        buffdata(1, :) = stim.yProbe;
        buffdata(2, :) = stim.ySupp;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        if k > stim.ThrowAway
            BothBuffs(k - stim.ThrowAway,  :) = vins;
        end
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
    fname = strcat(respDir,'SFOAE_linear_',subj,earname,'_',datetag, '.mat');
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

