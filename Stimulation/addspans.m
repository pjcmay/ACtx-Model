function [INP, linp] = addspans(INP,Break)
% Defines breaks in the simulation at the onset of each stimulus
ntype = length(INP);
linp = zeros(ntype,1);
for itype = 1:ntype
    [~,inpspan] = size(INP(itype).inp);
    INP(itype).inpspan = inpspan;
    if isstruct(Break)
        INP(itype).break = Break(itype).b;
    else
        INP(itype).break = Break;
    end
    linp(itype) = inpspan;
end
