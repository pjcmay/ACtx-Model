%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This allows one to regenrate the results shown in Figure 2 of May (2021).
% The model (Mod), and simulation parameters (SimCond) are stored in Fig2x
% files in the folder Data. After loading Fig2x, you can rerun the
% simulation experiment using the Exp_xxxxx functions.
% Patrick J. C. May, Lancaster University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentpath = cd;
addpath(genpath(currentpath));

%% Fig. 2A: The omission MMN
load Fig2A
Ytest = Exp_OmissionMMN(Mod, SimCond);


%% Fig. 2A: The stimulus repetition MMN
load Fig2B
Exp_RepetitionMMN(Mod, SimCond);


%% Fig. 2C: The global deviance MMN
load Fig2C
Exp_GlobalDev(Mod, SimCond);


%% Fig. 2D: The MMN & the multistandard control condition
load Fig2D
Exp_MultiStdControl(Mod, SimCond);


%% Fig. 2E: The MMN for deviants in an anisochronous sequence
load Fig2E
Exp_Anisochron(Mod, SimCond);
