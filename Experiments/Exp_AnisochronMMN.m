function Res = Exp_AnisochronMMN(Mod, SimCond, Res)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define stimuli - from Schwartze et al. (2012)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stimparam.dur = 300; % tone duration (ms) 
Stimparam.onoff = 5; % onset and offset ramp (ms)
Stimparam.basel = 50; % baseline
Stimparam.minisi = 200; % minimum ISI
Stimparam.maxisi = 1000; % maximum ISI
Stimparam.pdev = 0.2; % deviant probability
Stimparam.n = 400; % number of tones in block
Stimparam.fstd = 6; % standard frequency
Stimparam.fdev = 9; % deviant frequency
if nargin == 2
    nUnits = Mod.nUnits;
    ctxdelay = Mod.ctxdelay;
    dur = Stimparam.dur; % tone duration
    onoff = Stimparam.onoff; % onset and offset ramp
    STIM = makeTonestims(Mod.amp, dur, onoff, Mod.sigmain, Mod); % this makes 16 tone stims
    basel = Stimparam.basel;
    minisi = Stimparam.minisi-basel-ctxdelay;
    maxisi =  Stimparam.maxisi-basel-ctxdelay;
    postw = minisi-basel-ctxdelay;
    fstd = Stimparam.fstd;
    fdev = Stimparam.fdev;
    clear INP Break
    INP(802).inp = [];
    Break(802).b = [];
    INP(1).inp = [zeros(nUnits,basel+ctxdelay)...
        STIM(fstd).snd zeros(nUnits,postw)]; %#ok<*SAGROW>
    Break(1).b = basel+ctxdelay;
    INP(2).inp = [zeros(nUnits,basel+ctxdelay)...
        STIM(fdev).snd zeros(nUnits,postw)]; %#ok<*SAGROW>
    Break(2).b = basel+ctxdelay;
    if maxisi == minisi
        INP(3).inp = zeros(nUnits,minisi); %#ok<*SAGROW>
    else
        for i = 1:(maxisi-minisi)
            INP(2+i).inp = zeros(nUnits,i); %#ok<*SAGROW>
        end
    end
    [INP, ~] = addspans(INP,Break);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define sequences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SeqOB = makeSequence(Stimparam.n,Stimparam.pdev);
    Seq0 = [];
    for i = 1:length(SeqOB)
        if maxisi == minisi
            Seq0 = [Seq0 SeqOB(i) 3]; %#ok<*AGROW>
        else
            %y = 2;
            y = 2+randsample(maxisi-minisi-2,1);
            Seq0 = [Seq0 SeqOB(i) y+2];
        end
    end
    if isfield(SimCond,'Seq0')
        Seq0 = SimCond.Seq0;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run Simulations    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = simgate(INP, Seq0, Mod, SimCond);
    Res.Ysort = sortResponses(R,INP,Seq0,Mod);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gather Parameters and Results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Res.Stimparam = Stimparam;
    Res.INP = INP;
    Res.Seq0 = Seq0;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ysort = Res.Ysort;
Mod.MEGcompW = [-2 1 1 -2 1 0 0];
meg = sigmaMEG(Ysort, Mod.MEGcompW);

MEG1st = meg(1).megst;
MEG2st = meg(2).megst;
MEG1 = mean(MEG1st(:,2:end),2);
MEG2 = mean(MEG2st(:,1:end),2);
MEG1 = MEG1-mean(MEG1(1:Stimparam.basel));
MEG2 = MEG2-mean(MEG2(1:Stimparam.basel));

t = (1:length(MEG1))-Stimparam.basel;
figure(3); clf
ax = gca;
ax.FontSize = 18; 
hold on
plot(t,MEG1,'LineWidth',2);
plot(t,MEG2,'Color','r','LineWidth',2);
hold off
box on
legend('Standard','Deviant','FontSize', 20)
xlabel('Time (ms)','FontSize', 24)
ylabel('MEG (a.u.)','FontSize', 24)

axis([-50 350 -20 200])
