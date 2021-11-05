function Y = calcMEGcomp0(Y,Mod,INP,Seq0)
%% Calculates the MEG components on the basis of Y and deposits it back
% into Y. Weightings Mod.MEGw express the connection type and are 
% (arbitrarily) defined in makeWtypeXX.m. Component N = length(Mod.Wtype)+1.
% Note: last term is due to synaptic input to the system.

nTypes = length(Mod.Wtype);
MEG(nTypes+1).meg = [];
for ii = 1:length(Seq0)
    Faff = INP(Seq0(ii)).inp;
    Aaff = Y(ii).aaff;
    Syninp = zeros(size(Faff));
    temp = repmat(Aaff',Mod.nInArea,1);
    Syninp(1:Mod.nInArea*Mod.nIn,:) = temp;
    Syninp = Mod.waff*Syninp.*Faff;    
    U = Y(ii).u; % = U(:,1:nUnits);
    A = Y(ii).a; % = cortical adaptation;
    if isfield(Y,'ai')
        Ai = Y(ii).ai;
    else
        Ai = ones(size(A));
    end
    Ui = Y(ii).ui;
    for iConType = 1:nTypes-1
        MEG(iConType).meg = sum(Mod.Wtype(iConType).w .* ...
                                 Mod.Wee*(A.*gain(U,Mod.theta))',1)';
        for iField = 1:Mod.nAreas
            MEG(iConType).Field(iField).meg = sum(...
                              Mod.Field(iField).w.* ...
                              Mod.Wtype(iConType).w .* ...
                              Mod.Wee*(A.*gain(U,Mod.theta))'...
                                 ,1)';
        end
    end
    MEG(nTypes).meg = sum(Mod.Wtype(nTypes).w .*Mod.Wei*(Ai.*gain(Ui,Mod.theta))',1)';
    for iField = 1:Mod.nAreas
        MEG(nTypes).Field(iField).meg = sum(...
                              Mod.Field(iField).w.* ...
                              Mod.Wtype(nTypes).w .* ...
                              Mod.Wei*(Ai.*gain(Ui,Mod.theta))'...
                                 ,1)';
    end
    MEG(nTypes+1).meg = sum(Syninp)';
    %MEG = sum(WeeMEG*(A.*gain(U))',1);
    %MEGei = sum(Mod.megw(7)*Mod.Wei*(Ai.*gain(Ui))',1);
    %MEG = MEG + Mod.megw(6)*sum(Syninp) + MEGei;
    %Y(ii).meg = MEG';
    Y(ii).megcomp = MEG;
end