%% Analyze discrete sfoae

% load in stim result

%1000 hz
%load('/Volumes/USB DISK/SFOAE_discrete_supp_flip/Results/SFOAE_discrete_SHLEar_25-Apr-2022_18_00_22.mat')

%4000 Hz
%load('/Volumes/USB DISK/SFOAE_discrete_supp_flip/Results/SFOAE_discrete_SHLEar_25-Apr-2022_18_12_17.mat')

% Subtract for sfoae
Vsfe = stim.ProbeBuffs + stim.SuppBuffs - stim.BothBuffs; 
 
% set matricies for final values
freqs = stim.freq_probe; 
magSFOAE_pmtm = zeros(1,length(freqs)); 
magNOISE_pmtm = zeros(1,length(freqs)); 
magSFOAE_fft = zeros(1,length(freqs)); 
magNOISE_fft = zeros(1,length(freqs)); 
phi = zeros(1,length(freqs));

% Artifact rejection and set mean SFOAE and NOISE
for l = 1:stim.points
    
    trials = stim.trials; 
    first = (l-1)*trials+1; 
    last = l*trials;
    f0 = stim.freq_probe(l);
    stimend = floor(stim.fs*stim.toneDur); 
    
    SFOAEtrials = Vsfe(first:last,1:stimend); 
    
    % Artifact Rejection
    energy = squeeze(sum(SFOAEtrials.^2, 2));
    good = energy < median(energy) + 2*mad(energy);
    SFOAEclean = mean(SFOAEtrials(good, :), 1);

    good_2x = good(1: 2*floor(numel(good) / 2));
    Noiseclean = mean(SFOAEtrials(good_2x(1:2:end), :), 1) - mean(SFOAEtrials(good_2x(2:2:end), :), 1);

    SFOAE = SFOAEclean * stim.VoltageToPascal * stim.PascalToLinearSPL; 
    NOISE = Noiseclean * stim.VoltageToPascal * stim.PascalToLinearSPL; 
 
    nfft = 2^ceil(log2(length(SFOAE))); 
    
    % pmtm
    [Pxx, fp] = pmtm(SFOAE, 1.5, nfft, stim.fs); 
    [Pxx_n, f_n] = pmtm(NOISE, 1.5, nfft, stim.fs); 
    
    [~,ind] = min(abs(fp-f0));  
    magSFOAE_pmtm(1,l) = pow2db(Pxx(ind)); 
    magNOISE_pmtm(1,l) = pow2db(Pxx_n(ind));
   
    % fft 
    fft_sfe = fft(SFOAE, nfft); 
    fft_n = fft(NOISE, nfft); 
    f = (0:(nfft-1))*stim.fs/nfft; 
    [~,ind_f0] = min(abs(f-f0)); 
    %plot all these ffts
    %figure; 
    %plot(f(ind_f0-100:ind_f0+200),abs(fft_sfe(ind_f0-100:ind_f0+200)), f(ind_f0), abs(fft_sfe(ind_f0)), 'b*', freqs(l),abs(fft_sfe(ind_f0)),'gs', [freqs(l)-47,freqs(l)-47], [0,abs(fft_sfe(ind_f0))]);
    
    magSFOAE_fft(1,l) = db(abs(fft_sfe(ind_f0))); 
    magNOISE_fft(1,l) = db(abs(fft_n(ind_f0))); 
    phi(1,l) = angle(fft_sfe(ind_f0));  
    
    % inner product
    t_model =  (0:(numel(SFOAE)-1))/stim.fs;
    w = hanning(numel(SFOAE))';
    wsn = w.*sin(2*pi*f0*t_model);
    wcs = w.*cos(2*pi*f0*t_model);
    magSFOAE_ip(l) = db(sqrt( sum(wsn .* SFOAE)^2 + sum(wcs .* SFOAE)^2));
    magNOISE_ip(l) = db(sqrt( sum(wsn .* NOISE)^2 + sum(wcs .* NOISE)^2));
    phi_ip(l) = -angle(1j.* sum(wsn .* SFOAE) + sum(wcs .* SFOAE) );
end

phase_fft = unwrap(phi)./(2*pi);
[P_fft, S_fft] = polyfit(freqs,phase_fft, 1); 
P_fft(1)

phase_ip = unwrap(phi_ip)./(2*pi); 
[P_ip, S_ip] = polyfit(freqs,phase_ip, 1); 
P_ip(1)


% figure; 
% plot(freqs, magSFOAE_pmtm, 's-b', freqs, magNOISE_pmtm, 's-r');
% legend('SFOAE', 'NOISE')
% title('PMTM'); 
% figure; 
% subplot(2,1,1); 
% plot(freqs,magSFOAE_fft, 's-b', freqs, magNOISE_fft, 's-r'); 
% legend('SFOAE', 'NOISE')
% title('FFT');
% subplot(2,1,2); 
% plot(freqs,phase_fft); 
% title('phase fft'); 
figure; 
subplot(2, 1, 1); 
plot(freqs,magSFOAE_ip, 's-b', freqs, magNOISE_ip, 's-r'); 
legend('SFOAE', 'NOISE')
title('IP');
subplot(2, 1, 2); 
plot(freqs,phase_ip)
title('phase ip'); 
