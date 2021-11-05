function R = Exp_MultiStdControl(Mod, SimCond, R)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define stimuli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2
    nUnits = Mod.nUnits; % N columns/units
    ctxdelay = Mod.ctxdelay;
    dur = 50;  % tone duration (ms)
    onoff = 5; % onset and offset ramp (ms)
    STIM = makeTonestims(Mod.amp, dur, onoff, Mod.sigmain, Mod); % this makes 16 tone stims
    basel = 100; % baseline
    soa = 500; % stimulus-onset asynchrony
    isi = soa-basel-ctxdelay-dur; % offset-to-onset isi
    clear INP Break
    f = [10:-1:4 11:13]; % tone frequencies (correspondence with Nieto-Diego et al.)
    nfreq = length(f);
    INP(nfreq).inp = [];
    Break(nfreq).b = [];
    for i = 1:nfreq
        INP(i).inp = [zeros(nUnits,basel+ctxdelay) STIM(f(i)).snd zeros(nUnits,isi)]; %#ok<*SAGROW>
        Break(i).b = basel+ctxdelay;
    end
    [INP, ~] = addspans(INP,Break);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define sequences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nblock = SimCond.nblock;
    Seq = randperm(nfreq);
    if Seq(1) == 1  % assure that first stim is not a deviant
        Seq = fliplr(Seq);
    end
    for i = 1:nblock-1
        block = randperm(nfreq);
        if block(1) == Seq(length(Seq))  % make sure dev dev does not occur
            block = fliplr(block);
        end
        Seq = [Seq block]; %#ok<AGROW>
    end
    Seq0 = Seq;
    SeqControl = Seq0;
    SeqOddball = 1+(Seq0>1); % NOTE: Std = 2, Dev = 1
    SEQ(1).seq = SeqOddball; % Oddball sequence
    SEQ(2).seq = SeqControl; % Multistandard control sequence
    if isfield(SimCond,'SEQ')
        SEQ = SimCond.SEQ;
    end
    INPob = INP;
    INPctrl = INP;
    I(1).INP = INPob;
    I(2).INP = INPctrl;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run simulations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y(2).Y = [];
    tic
    for i = 1:2
        Y(i).Y = simgate(I(i).INP, SEQ(i).seq, Mod, SimCond);
    end
    toc
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sort responses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Yob = Y(1).Y;
    Ycnt = Y(2).Y;
    SeqOddball = SEQ(1).seq;
    SeqControl = SEQ(2).seq;
    Seq = SeqOddball;
    Seq(1:4) = 0; % Ignore initial responses
    Ysort_ob = sortResponses(Yob,INPob,Seq,Mod);
    Seq = SeqControl;
    Seq(1:4) = 0; % Ignore initial responses
    Seq(Seq>=2) = 2;
    Ysort_cnt = sortResponses(Ycnt,INPctrl,Seq,Mod);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gather parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Par.dur = dur;
    Par.onoff = onoff;
    Par.basel = basel;
    Par.soa = soa;
    Par.isi = isi;
    Par.f = f;
    Par.I = I;
    Par.nblock = nblock;
    Par.pdev = sum(SeqOddball == 1)/length(SeqOddball);
    R.Par = Par;
    R.SEQ = SEQ;
    R.Ysort_ob = Ysort_ob;
    R.Ysort_cnt = Ysort_cnt;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1
    MEGob = sigmaMEG(R.Ysort_ob, Mod.MEGcompW);
    MEGcnt = sigmaMEG(R.Ysort_cnt, Mod.MEGcompW);
    Mstd = MEGob(2).megm; % Note the switch from the norm: STIM(2) is the std.
    Mdev = MEGob(1).megm;
    Mms = MEGcnt(1).megm;
    Rstd = Mstd-mean(Mstd(1:R.Par.basel));
    Rdev = Mdev-mean(Mdev(1:R.Par.basel));
    Rms = Mms-mean(Mms(1:R.Par.basel));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEG Responses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1); clf
    ax = gca;
    ax.FontSize = 18; 
    lw1 = 2;
    taxis = linspace(-99,400,500);
    hold on
    plot(taxis,Rstd,'Color','b','LineWidth',lw1)
    plot(taxis,Rdev,'Color','r','LineWidth',lw1)
    plot(taxis,Rms,'Color','k','LineWidth',lw1)
    %plot(taxis,Rca,'Color','g','LineWidth',lw1)
    hold off
    box on
    legend('Standard','Deviant','Control','FontSize', 20)
    xlabel('Time (ms)','FontSize', 24)
    ylabel('MEG (a.u.)','FontSize', 24)
    axis([-50 350 -20 150])
    
 

end