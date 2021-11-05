function R = Exp_GlobalDev(Mod, SimCond, R)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define stimuli - from Wacongne et al. (2011)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2
    nUnits = Mod.nUnits;
    ctxdelay = Mod.ctxdelay;
    dur = 50; % tone duration (ms)
    onoff = 7; % onset and offset ramp (ms)
    STIM = makeTonestims(Mod.amp, dur, onoff, Mod.sigmain, Mod); % this makes 16 tone stims
    basel = 100; % baseline
    isi1 = 100;
    isi2 = 850-basel-ctxdelay;              % offset-to-onset isi
    f1 = 5;
    f2 = 12;
    n = 400;
    temp = zeros(size(STIM(1).snd));
    for i = 1:length(f1)
        temp = temp + STIM(f1(i)).snd;
    end
    S1 = temp;
    temp = zeros(size(STIM(1).snd));
    for i = 1:length(f1)
        temp = temp + STIM(f2(i)).snd;
    end
    S2 = temp;
    ISI1 = zeros(nUnits,isi1);
    ISI2 = zeros(nUnits,isi2);
    ISIexpected = zeros(nUnits,isi2+dur+isi1);
    npre = 10; % number of initial stimuli to discard
    pdev = 0.25; % deviant probability
    clear INP Break
    INP(2).inp = [];
    Break(2).b = [];
    
    % xxxx - sequence of 4 repeating tones
    INP(1).inp = [zeros(nUnits,basel+ctxdelay)...
        S1 ISI1 S1 ISI1 S1 ISI1 S1 ISIexpected]; %#ok<*SAGROW>
    Break(1).b = basel+ctxdelay + [0 dur+isi1 2*(dur+isi1) 3*(dur+isi1)];
    
    % xxxxY - sequence of 4 repeating tones followed by local deviant
    INP(2).inp = [zeros(nUnits,basel+ctxdelay)...
        S1 ISI1 S1 ISI1 S1 ISI1 S1 ISI1 S2 ISI2]; %#ok<*SAGROW>
    Break(2).b = basel+ctxdelay + [0 dur+isi1 2*(dur+isi1) 3*(dur+isi1)...
        4*(dur+isi1)];
    
    % xxxxX - sequence of 5 repeating tones
    INP(3).inp = [zeros(nUnits,basel+ctxdelay)...
        S1 ISI1 S1 ISI1 S1 ISI1 S1 ISI1 S1 ISI2]; %#ok<*SAGROW>
    Break(3).b = basel+ctxdelay + [0 dur+isi1 2*(dur+isi1) 3*(dur+isi1)...
        4*(dur+isi1)];
    
    [INP, ~] = addspans(INP,Break);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define sequences - from Wacongne et al. (2011)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Seq0 = makeSequence(n,pdev);
    Seq0 = [ones(1,npre) Seq0];
    devs1 = [3*ones(1,15) ones(1,10)];
    devs1 = repmat(devs1,1,n/100);
    devs1 = devs1(randperm(length(devs1)));
    devs2 = devs1;
    devs2(devs2==3) = 2;
    SeqOBxxxxY = [];
    SeqOBxxxxX = [];
    devindex = 1;
    for i = 1:length(Seq0)
        switch Seq0(i)
            case 1 % Define standard
                SeqOBxxxxY(i) = 2;  %#ok<*AGROW>
                SeqOBxxxxX(i) = 3;
            case 2 % Define deviants
                SeqOBxxxxY(i) = devs1(devindex);
                SeqOBxxxxX(i) = devs2(devindex);
                devindex = devindex+1;
        end
    end
    % Condition 1: xxxx occurs 100%
    SEQ(1).seq = ones(1,length(Seq0));
    % Condition 2: xxxxY [75%], xxxxX [15%], xxxx [10%]
    SEQ(2).seq = SeqOBxxxxY;
    % Condition 3: xxxxY [15%], xxxxX [75%], xxxx [10%]
    SEQ(3).seq = SeqOBxxxxX;
    % Note: percentages above pertain to situation where the initial npre
    % tones have been discarded
    if isfield(SimCond,'SEQ')
        SEQ = SimCond.SEQ;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run simulations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y(3).Y = [];
    tic
    for i = 1:3
        Y(i).Y = simgate(INP, SEQ(i).seq, Mod, SimCond);
    end
    toc
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort responses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    YExp = Y(1).Y; % Responses in Condition 1 - xxxx 100%
    YxxY = Y(2).Y; % Responses in Condition 2 - xxxxY 75%
    YxxX = Y(3).Y; % Responses in Condition 3 - xxxxX 75%
    SeqExp = SEQ(1).seq;
    SeqxxY = SEQ(2).seq;
    SeqxxX = SEQ(3).seq;
    Seq = SeqExp;
    Seq(1:npre) = 0; % discard npre initial responses
    Ysort_exp = sortResponses(YExp,INP,Seq,Mod); % Condition 1
    Seq = SeqxxY;
    Seq(1:npre) = 0;
    Ysort_xxY = sortResponses(YxxY,INP,Seq,Mod); % Condition 2
    Seq = SeqxxX;
    Seq(1:npre) = 0;
    Ysort_xxX = sortResponses(YxxX,INP,Seq,Mod); % % Condition 3
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gather parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Par.dur = dur;
    Par.onoff = onoff;
    Par.basel = basel;
    Par.isi1 = isi1;
    Par.isi2 = isi2;
    Par.f1 = f1;
    Par.f2 = f2;
    Par.n = n;
    Par.npre = npre;
    Par.pdev = pdev;
    R.Mod = Mod;
    R.Par = Par;
    R.SEQ = SEQ;
    R.Ysort_exp = Ysort_exp;
    R.Ysort_xxY = Ysort_xxY;
    R.Ysort_xxX = Ysort_xxX;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Get defintions
% -------------------------------------------------------------------------
basel = R.Par.basel;
dur = R.Par.dur;
isi1 = R.Par.isi1;
FLT = 1; % filter switch
FcHi = 0.1; % (Hz)
Fs    = 1000; % sampling rate

t0  = basel+[0 dur+isi1 2*(dur+isi1) 3*(dur+isi1)...
    4*(dur+isi1)]; % tone onsets
% -------------------------------------------------------------------------
% Get ERFs
% -------------------------------------------------------------------------
% Responses to xxxxX
Rstd = MEGxxX(3).megm; % xxxxX as standard
Rdev = MEGxxY(3).megm; % xxxxX as deviant (when xxxxY is standatd)

% Filtering and baseline correction
if FLT == 1
    Rstd = highpass(Rstd,FcHi,Fs);
    Rdev = highpass(Rdev,FcHi,Fs);
end
Rstd = Rstd-mean(Rstd(1:R.Par.basel));
Rdev = Rdev-mean(Rdev(1:R.Par.basel));
% -------------------------------------------------------------------------
% Plot ERFs
% -------------------------------------------------------------------------
t = (1:length(Rexp))-t0(5); % set zero time to onset of 5th tone
figure(1); clf
sc = 0.63;
ax = gca;
ax.FontSize = sc*18;
%title('Unexpected repetition')
hold on
plot(t,Rstd,'LineWidth',2);
plot(t,Rdev,'Color','r','LineWidth',2);
hold off
box on
pbaspect([4 1 1])
legend('Frequent','Rare','FontSize', sc*20)
xlabel('Time (ms)','FontSize', sc*24)
ylabel('MEG','FontSize', sc*24)
axis([-700 500 -20 250])
