% change these vars if needed
windowdur = 0.04; % 40ms in paper
offsetwin = 0.02; % 20ms in paper
npoints = 32; 

% range that can be tested
first = ceil((windowdur/2)*stim.Fs); % sample
last = numel(stim.t) - first - ceil(offsetwin*stim.Fs); % sample

% setting stuff from stim
phiProbe_inst = stim.phiProbe_inst; 
fmin = stim.fmin; 
fmax = stim.fmax; 
t = stim.t; 
    
% set SFOAE
SFOAEtrials = stim.ProbeBuffs + stim.SuppBuffs - stim.BothBuffs; 
% SFOAEtrials = stim.ProbeBuffs;
% Artifact rejection and set mean SFOAE and NOISE
energy = squeeze(sum(SFOAEtrials.^2, 2));
good = energy < median(energy) + 2*mad(energy);
SFOAEclean = mean(SFOAEtrials(good, :), 1);

good_2x = good(1: 2*floor(numel(good) / 2));
Noiseclean = mean(SFOAEtrials(good_2x(1:2:end), :), 1) - mean(SFOAEtrials(good_2x(2:2:end), :), 1);
    
SFOAE = SFOAEclean; % * stim.VoltageToPascal * stim.PascalToLinearSPL; 
NOISE = Noiseclean; % * stim.VoltageToPascal * stim.PascalToLinearSPL;

if stim.speed > 0 % upsweep
    freq_inst = fmin + abs(stim.speed)*t; 
    testfreq = linspace(freq_inst(1,first), freq_inst(1,last), npoints); 
    t_freq = (testfreq-f1)/stim.speed; 

else % downsweep
    freq_inst = fmax - abs(stim.speed)*t; 
    testfreq = linspace(freq_inst(1,first), freq_inst(1,last), npoints); 
    t_freq = (fmax-testfreq)/abs(stim.speed); 
end

% % set freq we're testing and the timepoints when they happen.  
% freq_inst = f1 + stim.speed*t; 
% %npoints = ceil(abs((freq_inst(first) - freq_inst(last))) / (windowdur* abs(stim.speed))); 
% testfreq = linspace(freq_inst(1,first), freq_inst(1,last), npoints); 
% t_freq = (testfreq-f1)/stim.speed; 

% Set empty matricies for next steps
maxoffset = ceil(stim.Fs * offsetwin);
coeffs = zeros(npoints, 2); 
coeffs_n = zeros(npoints, 2); 
tbest = zeros(npoints, 1); 
tbest_n = zeros(npoints, 1); 

% Generate model of chirp and test against response
for k = 1:npoints 
    fprintf(1, 'Running window %d / %d\n', k, (npoints));
    
    win = find( (t > (t_freq(k)-windowdur/2)) & (t < t_freq(k)+windowdur/2)); 
    taper = hanning(numel(win))';
       
    model = [cos(2*pi*phiProbe_inst(win)) .* taper; 
        -sin(2*pi*phiProbe_inst(win)) .* taper]; 
   
    
    % zero out variables for offset calc
    coeff = zeros(maxoffset, 2); 
    coeff_n = zeros(maxoffset, 2); 
    resid = zeros(maxoffset, 1); 
    resid_n = zeros(maxoffset, 1); 
    
     for offset = 0:maxoffset
        resp = SFOAE(win+offset); % .* taper; 
        resp_n = NOISE(win+offset); % .* taper;
         
        coeff(offset + 1, :) = model' \ resp'; 
        coeff_n(offset + 1, :) = model' \ resp_n';  
         
        resid(offset + 1) = sum( (resp - (coeff(offset+1, :) * model)).^2); 
            
     end
     
  
    [~, ind] = min(resid);  

     coeffs(k, :) = coeff(ind, :); 
     coeffs_n(k, :) = coeff_n(ind, :); 
     
     tbest(k) = (ind - 1) * (1/stim.Fs); % delay in sec
    
     
end

a = coeffs(:, 1); 
b = coeffs(:, 2); 
a_n = coeffs_n(:, 1); 
b_n = coeffs_n(:, 2); 

oae = abs(complex(a, b)) .* stim.VoltageToPascal .* stim.PascalToLinearSPL; 
nf = abs(complex(a_n, b_n)) .* stim.VoltageToPascal .* stim.PascalToLinearSPL;  

phi = -tbest .* testfreq' ;  % + angle(complex(a,b))/(2*pi); 


% Plot figures
figure;
subplot(2,1,1)
plot(testfreq, db(oae), '-s',  testfreq, db(nf), '--', 'linew', 2); 
legend('OAE', 'Noise', 'Location', 'southeast'); 
ylabel('SFOAE level (dB SPL)', 'FontSize', 16);

subplot(2,1,2)
plot(testfreq, unwrap(phi), '-o');
xlabel('Probe Frequency (Hz)', 'FontSize', 16); 
ylabel('Phase (cycles)', 'FontSize', 16); 