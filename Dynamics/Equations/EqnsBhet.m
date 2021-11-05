function dy = EqnsBhet(t, y, Mod, Inputs)
% These equations include heterogeneous adaptation of the Wie synapses
Syninput = Inputs.Syninp;
%Fpfc = Inputs.pfc;
nUnits = Mod.nUnits;
Wee = Mod.Wee;
Wei = Mod.Wei;
Wie = Mod.Wie;
Wii = Mod.Wii;
tau = Mod.tau; % membrane constant, excitatory [N]
tau_i = Mod.tau_i; % membrane constant, inhibitory [N]
theta = Mod.theta;

set1 = Mod.adapset1; % adaptation constants in MGV to be used for this set of
% cortical columns -> indices of columns stored in row vector Mod.adapset1 [N]

set2 = Mod.adapset2; % adaptation constants in Cortex to be used for this set of
% cortical columns -> indices of columns stored in row vector Mod.adapset2 [N]

taua1 = Mod.taua1; % adaptation decay time constant in MGV [N]
tauaon1 = Mod.tauaon1; % adaptation onset time constant in MGV [N]

taua2 = Mod.taua2; % adaptation decay time constant in Cortex [N]
tauaon2 = Mod.tauaon2; % adaptation onset time constant in Cortex [N]

ue = y(1:nUnits); % range of y for excitatory state variable[N]
ui = y(nUnits+1:2*nUnits); % range of y for inhibitory state variable [N]
a1 = y(2*nUnits+1:3*nUnits); % range of y for adaptation in MGV [N]
a2 = y(3*nUnits+1:4*nUnits); % range of y for adaptation in Cortex [N]


% Excitation
% in MGV [N]
dy(set1) = 1/tau * (-ue(set1) + (Wee(set1,:)*(a1.*gain(ue,theta)))  + ...
                           -Wei(set1,:)*gain(ui,theta) + ...
                           Mod.waff*Syninput(set1,int64(t*1000)) + ...
                           Mod.tonic(set1));
% Inhibition
% in MGV [N]
dy(nUnits+set1) = (1/tau_i)* (-ui(set1)) + ...
                  (1/Mod.tau_i_on)*(Wie(set1,:)*(a1.*gain(ue,theta)) +...
                                   -Wii(set1,:)*gain(ui,theta));


% Excitation
% in Cortex [N]
dy(set2) = 1/tau * (-ue(set2) + (Wee(set2,:)*(a2.*gain(ue,theta)))  + ...
                           -Wei(set2,:)*gain(ui,theta) + ...
                           Mod.waff*Syninput(set2,int64(t*1000)) ...
                          + Mod.tonic(set2));
% Inhibition
% in Cortex [N]
dy(nUnits+set2) = (1/tau_i)* (-ui(set2)) + ...
                  (1/Mod.tau_i_on)*(Wie(set2,:)*(a2.*gain(ue,theta)) + ...
                                   -Wii(set2,:)*gain(ui,theta));
                   
                       
                       
% Adaptation

% in MGV
dy(2*nUnits+1:3*nUnits) = (1 - a1)/taua1 ...
                           - (a1.*gain(ue,theta))/tauaon1;
% in Cortex                       
dy(3*nUnits+1:4*nUnits) = (1 - a2)/taua2 ...
                           - (a2.*gain(ue,theta))/tauaon2;
                       
                       
dy = dy';

