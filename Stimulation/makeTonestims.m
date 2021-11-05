function STIM = makeTonestims(amp, dur, onoff, sigma, Mod)
%   Makes "pure tone" stimuli for model
%   dur - duration of tones
%   onoff - duration of onset and offset ramps
%   sigma - width of tuning
%   nTargA - number of (core) areas targeted by input
%   Copyright Patrick May, Lancaster University, 2021

% Define parameters:

nUnits = Mod.nUnits;  % total number of units/columns in model

nIn = Mod.nIn; % number of columns per input area (16)

% Creating on and off ramp [N]

a = ones(1,dur); % Row vector of length 'duration', filled with ones [N]

a(1:onoff+1) = linspace(0,1,onoff+1); % content of entries 1 to onoff+1 are
% evenly spaced values between 0 and 1 [N]

a(dur-onoff:end) = linspace(1,0,onoff+1); % content of entries (dur-onoff)
% until end of vector (i.e. entry at index dur) are evenly spaced values 
% between 1 and 0 [N]


STIM(nIn).snd = [];

for f = 1:nIn
    
    profile = amp*normpdf(1:nIn,f,sigma)/max(normpdf(1:nIn,f,sigma));
    % normal distribution centred on column f [N]
    
    profile(profile<0.01) = 0; % zero input to columns outside tuning
    
    temp = repmat(a,nIn,1); % repeat copies of row vector 'a' to form
    % nIn x length(a) matrix [N]
    
    profile = repmat(profile',1,dur); % repeat copies of column vector
    % profile' to form a length(profile) x dur matrix [N]
    
    snd = temp.*profile; % element-wise multiplication of temp and profile,
    % strength of input modulated by tuning curve [N]
    
    temp = zeros(nUnits,dur);
    
    for iIn = 1:length(Mod.InArea) % Loop through all input areas [N]
        temp(Mod.Area(Mod.InArea(iIn)).icol,1:dur) = Mod.InModif(iIn)*snd;
    end
    STIM(f).snd = temp;
end

%% Random, multipeaked input
% Mod.InTonot: if 1, then tonotopic, if > 1 then n = number of random freqs
% targeted per input location. [N]
[~,indRandAreas] = find(Mod.InTonot > 1);
STIM0 = STIM;
for iRA = 1:length(indRandAreas)
    RND = Mod.InRand(iRA).RND;
    [~,nfrand] = size(RND);
    for f = 1:nIn
        temp = zeros(size(STIM0(1).snd));
        for i = 1:nfrand
            frand = RND(f,i);
            temp = temp+STIM0(frand).snd;
        end
        ColumnsInQuestion = Mod.Area(Mod.InArea(indRandAreas(iRA))).icol;
        STIM(f).snd(ColumnsInQuestion,:) = temp(ColumnsInQuestion,:);
    end
end


