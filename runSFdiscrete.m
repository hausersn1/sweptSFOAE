try 
    % Initialize ER-10X  (Also needed for ER-10C for calibrator)
    initializeER10X;
    
    % Initializing TDT
    % Specify path to cardAPI here
    pcard = genpath('C:\Experiments\cardAPI\');
    addpath(pcard);
    card = initializeCard;
    
    
    % Get stimulus structure
    stim = makeSFdiscrete();
    
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
  %% Present stimuli 1 trial at a time 
  
  ProbeBuffs = zeros(stim.trials*stim.points, size(stim.t,2)); 
  SuppBuffs = zeros(stim.trials*stim.points, size(stim.t,2)); 
  BothBuffs = zeros(stim.trials*stim.points, size(stim.t,2)); 
  
  for k=1:stim.points 
      for j = 1:stim.trials
          
        delayComp = 1; % Always
        
        % Do probe only
        dropSupp = 120;
        dropProbe = stim.drop_probe;
        stims_probe = scaleSound(rampsound(cos(2*pi*stim.freq_probe(k)*stim.t), stim.fs, 0.005)); 
        buffdata = zeros(2, (numel(stim.t)+stim.sampWaitDur));
        buffdata(1, :) = stims_probe;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        
        ProbeBuffs((k-1)*stim.trials+j,:) = vins; 
        WaitSecs(max(stim.t)+0.05); 
        
        % Do suppressor only
        dropProbe = 120;
        dropSupp = stim.drop_supp;
        stims_supp = scaleSound(rampsound(cos(2*pi*stim.freq_supp(k)*stim.t+stim.phi(j)), stim.fs, 0.005)); 
        buffdata = zeros(2, (numel(stim.t)+stim.sampWaitDur));
        buffdata(2, :) = stims_supp;
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        
        SuppBuffs((k-1)*stim.trials+j,:) = vins; 
        WaitSecs(max(stim.t)+0.05); 
        
        % Do both 
        dropSupp = stim.drop_supp;
        dropProbe = stim.drop_probe;
        stims_probe = scaleSound(rampsound(cos(2*pi*stim.freq_probe(k)*stim.t), stim.fs, 0.005)); 
        stims_supp = scaleSound(rampsound(cos(2*pi*stim.freq_supp(k)*stim.t+stim.phi(j)), stim.fs, 0.005));
        buffdata = zeros(2, (numel(stim.t)+stim.sampWaitDur));
        buffdata(1, :) = stims_probe;
        buffdata(2, :) = stims_supp; 
        vins = playCapture2(buffdata, card, 1, 0,...
            dropProbe, dropSupp, delayComp);
        
        BothBuffs((k-1)*stim.trials+j,:) = vins; 
        WaitSecs(max(stim.t)+0.05); 
        

      end
      
      fprintf(1, 'Done with frequency %d / %d\n', k, stim.points)
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
  fname = strcat(respDir,'SFOAE_discrete_',subj,earname,'_',datetag, '.mat');
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
%% Stim structure (probebuffs, suppbuffs, both buffs)
% if 10 trials and three frequencies tested
% row 1-10 = each trial of first freq
% row 11-20 = each trial of second freq
% row 21-30 = each trial of third freq
  