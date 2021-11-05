function [R, Init] = simgate(INP, STIMSEQ, Mod, SimCond, Init)
% This picks the right equations to use in the simulations, depending on
% the value of Sim.mode. This function is redundant in this release because
% only the equations for Mod.mode = 5 are provided. These are the ones
% corresponding to those given in the paper May (2021).

if nargin == 4
    Init = [];end
switch Mod.mode
    case 5  
        [R, Init] = runsimB(INP, STIMSEQ, Mod, SimCond, Init);        
end
