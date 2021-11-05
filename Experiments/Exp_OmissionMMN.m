function Res = Exp_OmissionMMN(Mod, SimCond, Res)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stimparam.dur = 50; % tone duration
Stimparam.onoff = 5; % onset and offset ramp
Stimparam.basel = 0; % baseline
Stimparam.soa = 100; % From Yabe
Stimparam.pdev = 0.5; % Probability of deviant
Stimparam.n = 400; % Number of stimuli in sequence
Stimparam.n = 10; % Number of stimuli in sequence
Stimparam.f = 7; % Frequency of stimulus (1...16)
if nargin == 2
    nUnits = Mod.nUnits;
    dur = Stimparam.dur; % tone duration
    onoff = Stimparam.onoff; % onset and offset ramp
    STIM = makeTonestims(Mod.amp, dur, onoff, Mod.sigmain, Mod); % this makes 16 tone stims
    soa = Stimparam.soa;
    isi = soa-dur;              % offset-to-onset isi
    f = Stimparam.f;
    clear INP Break
    INP(2).inp = [];
    Break(2).b = [];
    % Construct sequence of 5 standards [S S S S S]
    INP(1).inp = [zeros(nUnits,isi) STIM(f).snd zeros(nUnits,isi)...
                                    STIM(f).snd zeros(nUnits,isi)...
                                    STIM(f).snd zeros(nUnits,isi)...
                                    STIM(f).snd zeros(nUnits,isi)...
                                    STIM(f).snd]; %#ok<*SAGROW>
    Break(1).b = [isi (2*isi+dur) (3*isi+2*dur) (4*isi+3*dur) (5*isi+4*dur)];
    % Construct sequence with omission [S Om S S S]
    INP(2).inp = [zeros(nUnits,isi) STIM(f).snd zeros(nUnits,isi)...
                                    0*STIM(f).snd zeros(nUnits,isi)...
                                    STIM(f).snd zeros(nUnits,isi)...
                                    STIM(f).snd zeros(nUnits,isi)...
                                    STIM(f).snd]; %#ok<*SAGROW>
    Break(2).b = [isi (2*isi+dur) (3*isi+2*dur) (4*isi+3*dur) (5*isi+4*dur)];
    [INP, ~] = addspans(INP,Break);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define sequences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(SimCond,'Seq0')
        Seq0 = SimCond.Seq0;
    else
        Seq0 = makeSequence(Stimparam.n,Stimparam.pdev);
    end
    display(sum(Seq0==2)/length(Seq0))
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run Simulations up to Permutation Point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = simgate(INP, Seq0, Mod, SimCond);
    %%
    Res.Ysort = sortResponses(R,INP,Seq0,Mod);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gather Parameters and Results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Res.Stimparam = Stimparam;
    Res.INP = INP;
    Res.Seq0 = Seq0;
    Res.SimCond = SimCond;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INP = Res.INP;
t0 = INP(1).break; % Stimulus onsets
Ysort = Res.Ysort;
Mod.MEGcompW = [-2 1 1 -2 1 0 0]; % Connection type weights for MEG
meg = sigmaMEG(Ysort, Mod.MEGcompW);

MEG1st = meg(1).megst; % single-trial responses std
MEG2st = meg(2).megst; % single-trial responses omission
MEG1 = mean(MEG1st(:,3:end),2); % ERF - discard intial responses
MEG2 = mean(MEG2st(:,1:end),2); % ERF containing omission response

% Construct ERFs in time window extending from -50 ms to 350 ms relative
% to stimulus onset. This can be done for the 1st and 2nd stimulus of the 
% sequence (because the total time window of sequence is 500 ms). This
% means that the ERF for the omission captures all the stimulus omissions
% and the ERF to the standards captures the first and second stimulus of
% the 5-stimulus sequence.
span1 = (t0(1)-49):(t0(1)+350); % [-50 350] ms relative to onset of stim 1
span2 = (t0(2)-49):(t0(2)+350); % [-50 350] ms relative to onset of stim 2
Std1 = MEG1(span1); % ERF to first stimulus in sequence
Std2 = MEG1(span2); % ERF to second stimulus in sequence
Std = mean([Std1 Std2],2); % Averaging across 1st and 2nd stimulus
Dev = MEG2(span2);

% Filtering and baseline correction
FcHi = 1; % Hz
Fs    = 1000;
Std = highpass(Std,FcHi,Fs);
Dev = highpass(Dev,FcHi,Fs);
Std = Std-mean(Std(1:100));
Dev = Dev-mean(Dev(1:100));

% Single-trial responses
figure(1); clf
subplot(2,1,1)
plot(MEG1st)
title('Single trial - standards')
subplot(2,1,2)
plot(MEG2st)
title('Single trial - omission')

% Fig 2A
t = (1:length(Std))-50;
figure(2); clf
ax = gca;
ax.FontSize = 18; 
hold on
plot(t,Std,'LineWidth',2);
plot(t,Dev,'Color','r','LineWidth',2);
hold off
box on
legend('Standard','Omission','FontSize', 20)
xlabel('Time (ms)','FontSize', 24)
ylabel('MEG','FontSize', 24)
axis([-50 350 -15 35])
