%% Stimulus Generation for Discrete SFOAE 
function stim = makeSFdiscrete() 
%% Set Variable options
fc = 4000; 
trials = 64; % trial = probe, supp, probe+supp
resol = 32; % Hz between test frequencies was 47
points = 10; % odd to be centered; was 17
toneDur = 0.125; % seconds per stim 
waitDur = 0.0; % seconds at end of stim was 20 ms
fs = 48828.125; % samples per second
drop_probe = 58; % attenuation of probe  53
drop_supp = 38; % attenuation of suppressor 33 
supp_addHz = -47; % Hz added to probe frequency

%% Create matrices of probe and supp tones

t = 0:(1/fs):(toneDur-1/fs);
t_wait = round(fs*waitDur); 

testwidth = (points-1)*resol; 
freq_probe = (fc - (testwidth/2)):resol:(fc + (testwidth/2)); 
freq_supp = freq_probe + supp_addHz; 


phi = 2*pi*(0:(trials-1))/trials; 

% for i = 1:length(freq_probe) 
%     stims_probe(i,1:length(t)) = scaleSound(rampsound(cos(2*pi*freq_probe(i)*t), fs, 0.005)); 
%     stims_supp(i,1:length(t)) = scaleSound(rampsound(cos(2*pi*freq_supp(i)*t+phi), fs, 0.005)); 
% end

%% Save all to stim struct
    
stim.fc = fc; 
stim.trials = trials; 
stim.points = points; 
stim.fs = fs; 
stim.t = t; 
stim.freq_probe = freq_probe; 
stim.freq_supp = freq_supp; 
%stim.stims_supp = stims_supp; stim.stims_probe = stims_probe; 
stim.phi = phi; 
stim.drop_probe = drop_probe; 
stim.drop_supp = drop_supp; 
stim.toneDur = toneDur; 
stim.waitDur = waitDur; 
stim.sampWaitDur = t_wait; 

