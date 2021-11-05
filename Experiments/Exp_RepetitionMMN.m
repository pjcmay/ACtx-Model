function Res = Exp_RepetitionMMN(Mod, SimCond, Res)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stimparam.dur = 50; % tone duration (ms)
Stimparam.onoff = 5; % onset and offset ramp
Stimparam.basel = 0; % baseline
Stimparam.soa = 500; % stimulus onset asynchrony (ms)
Stimparam.pdev = 0.2; % deviant probability
Stimparam.n = 200; % Number of stimulus sequences
Stimparam.f = [6 9]; % frequencies for A and B
if nargin == 2
    nUnits = Mod.nUnits;
    dur = Stimparam.dur; % tone duration
    onoff = Stimparam.onoff; % onset and offset ramp
    STIM = makeTonestims(Mod.amp, dur, onoff, Mod.sigmain, Mod); % this makes 16 tone stims
    soa = Stimparam.soa;
    isi = soa-dur;              % offset-to-onset isi
    f1 = Stimparam.f(1);
    f2 = Stimparam.f(2);
    clear INP Break
    INP(2).inp = [];
    Break(2).b = [];
    % Construct sequence ABAB
    INP(1).inp = [zeros(nUnits,isi/2) STIM(f1).snd zeros(nUnits,isi)...
                                      STIM(f2).snd zeros(nUnits,isi)...
                                      STIM(f1).snd zeros(nUnits,isi)...
                                      STIM(f2).snd zeros(nUnits,isi/2)]; %#ok<*SAGROW>
    Break(1).b = [isi/2 isi/2+dur+isi (isi/2+2*dur+2*isi) (isi/2+3*dur+3*isi)];
    % Construct sequence AAAB
    INP(2).inp = [zeros(nUnits,isi/2) STIM(f1).snd zeros(nUnits,isi)...
                                      STIM(f1).snd zeros(nUnits,isi)...
                                      STIM(f1).snd zeros(nUnits,isi)...
                                      STIM(f2).snd zeros(nUnits,isi/2)]; %#ok<*SAGROW>
    Break(2).b = [isi/2 isi/2+dur+isi (isi/2+2*dur+2*isi) (isi/2+3*dur+3*isi)];
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
    % Run Simulations
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
t0 = INP(1).break;
Ysort = Res.Ysort;
Mod.MEGcompW = [-2 1 1 -2 1 0 0]; % Connection type weights for MEG
meg = sigmaMEG(Ysort, Mod.MEGcompW);

MEG1st = meg(1).megst; % single-trial responses std
MEG2st = meg(2).megst; % single-trial responses dev

MEG1 = mean(MEG1st(:,3:end),2); % ERF - discard intial responses
MEG2 = mean(MEG2st(:,1:end),2); % ERF - containing repetition

span1 = (t0(1)-99):(t0(1)+500); % time window for tone 1
span2 = (t0(2)-99):(t0(2)+500); % time window for tone 2 (repetition)
span3 = (t0(3)-99):(t0(3)+500); % time window for tone 3
Std1 = MEG1(span1);
Std2 = MEG1(span3);
Std = mean([Std1 Std2],2); % ERF to tone A in std sequence ABAB
Dev = MEG2(span2); % ERF to tone A when an unexpected repetition

% Filtering and baseline correction
FcHi = 1;
Fs    = 1000;
Std = highpass(Std,FcHi,Fs);
Dev = highpass(Dev,FcHi,Fs);
Std = Std-mean(Std(1:100));
Dev = Dev-mean(Dev(1:100));

% Plot 2B
t = (1:length(Std))-100;
figure(1); clf
ax = gca;
ax.FontSize = 18; 
hold on
plot(t,Std,'LineWidth',2);
plot(t,Dev,'Color','r','LineWidth',2);
hold off
box on
legend('Standard','Repetition','FontSize', 20)
xlabel('Time (ms)','FontSize', 24)
ylabel('MEG','FontSize', 24)
axis([-50 350 -20 110])


