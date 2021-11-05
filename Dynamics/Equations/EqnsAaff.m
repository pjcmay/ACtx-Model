function dy = EqnsAaff(t, y, Mod, inpaff)
% These equations determine afferent synaptic strengths

%dy = (1 - y)/Mod.tau_adap_r - Mod.a_adap*(y.*inpaff(:,int64(t*1000)));

dy = (1 - y)/Mod.tau_aff - (y.*inpaff(:,int64(t*1000)))/Mod.tauon_aff;



