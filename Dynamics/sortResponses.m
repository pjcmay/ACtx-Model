function Ysort = sortResponses(Y,INP,Seq,Mod,Sortparam)
%% Sorts responses according to stimulus identifier in Seq
% No baseline corrections or filtering here.
% Patrick J. C. May
if nargin == 5
    maxwin = Sortparam.maxwin;
else
    maxwin = inf;
end
ntype = max(Seq);
if isfield(Y,'megcomp')
    ncomp = length(Y(1).megcomp);
if isfield(Y(1).megcomp(1),'Field')
    nField = length(Y(1).megcomp(1).Field);
end
end
Ysort(ntype).placeholder = [];
for itype = 1:ntype;
    MEGc = [];
    win = min([INP(itype).inpspan maxwin]);
    U = zeros(win, Mod.nUnits, sum(Seq == itype));
    A = zeros(win, Mod.nUnits, sum(Seq == itype));
    Ui = zeros(win, Mod.nUnits, sum(Seq == itype));
    bingo = find(Seq == itype);
    for j = 1:sum(Seq == itype)
        if isfield(Y,'megcomp')
            for icomp = 1:ncomp
                temp = Y(bingo(j)).megcomp(icomp).meg;
                temp = temp(1:win)';
                MEGc(icomp).meg(:,j) = temp;  %#ok<AGROW>
                if isfield(Y(1).megcomp(1),'Field') && icomp<ncomp %FIX!!!! calcMEGcomp0 - afferent input not included now
                    tempMEG = [];
                    for iField = 1:nField
                        temp = Y(bingo(j)).megcomp(icomp).Field(iField).meg;
                        temp = temp(1:win)';
                        tempMEG = [tempMEG; temp]; %#ok<AGROW>
                    end
                    MEGc(icomp).megf(:,:,j) = tempMEG; %#ok<AGROW>
                end
            end
        end
        if isfield(Y,'u')
            temp = Y(bingo(j)).u;
            temp = temp(1:win,:);
            U(:,:,j) = temp;
        end
        if isfield(Y,'a')
            temp = Y(bingo(j)).a;
            temp = temp(1:win,:);
            A(:,:,j) = temp;
        end
        if isfield(Y,'ui')
            temp = Y(bingo(j)).ui;
            temp = temp(1:win,:);
            Ui(:,:,j) = temp;
        end
    end    
    if isfield(Y,'megcomp')
        Ysort(itype).megcomp = MEGc;
    end
    if isfield(Y,'u')
        Ysort(itype).u = U;
        %Ysort(itype).um = mean(U,3);
    end
    if isfield(Y,'a')
        Ysort(itype).a = A;
        %Ysort(itype).am = mean(A,3);
    end
    if isfield(Y,'ui')
        Ysort(itype).ui = Ui;
        %Ysort(itype).uim = mean(Ui,3);
    end
end
end

