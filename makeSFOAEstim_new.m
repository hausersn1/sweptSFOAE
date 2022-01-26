function stim = makeSFOAEstim_new(rawstim)

% rawstim (structure) should contain fields fmin, fmax, speed, Fs, ratio,
% VtoPforH

if ~exist('rawstim', 'var')
    stim.cf = 1000;
    stim.fmin = stim.cf/sqrt(2); % 1/2 octave below
    stim.speed = -2000; %Hz per second
    stim.diff = 50; % Hz (Fprobe - 50 = Fsupp; Probe is higher)
    stim.Fs = 48828.125;
    stim.fmax = stim.cf*sqrt(2); % 1/2 octave above
else
    stim = rawstim;
end
stim.drop_Probe = 50;
stim.drop_Supp = 30;
stim.ThrowAway = 2;
stim.Averages = 128;

fmin = stim.fmin;
fmax = stim.fmax;
dur = abs(fmax - fmin) / abs(stim.speed);
Fs = stim.Fs;
t = 0: (1/Fs): (dur - 1/Fs);
stim.t = t;

if stim.speed < 0
    phiProbe_inst = fmax*t + stim.speed*t.^2/2; % Cycles
else
    phiProbe_inst = fmin*t + stim.speed*t.^2/2; % Cycles
end

phiSupp_inst = phiProbe_inst - stim.diff*t;
stim.yProbe = scaleSound(rampsound(cos(2 * pi * phiProbe_inst), stim.Fs, 0.005));
stim.ySupp = scaleSound(rampsound(cos(2 * pi * phiSupp_inst), stim.Fs, 0.005));
stim.phiProbe_inst = phiProbe_inst;
stim.phiSupp_inst = phiSupp_inst;
