function MEG = sigmaMEG(Ysort, MEGcompW)
%%  Calculate the total MEG on the basis of components
%   (1) Build single-trial responses.
%   (2) Build averages of each component.
%   Sum components according to Mod.MEGcompW defined outside.
%   Patrick J. C. May, Lancaster University
ncomp = length(MEGcompW);
nsort = length(Ysort);
MEG(nsort).megm = [];
for is = 1:nsort
    if ~isempty(Ysort(is).megcomp)
        [win, nstim] = size(Ysort(is).megcomp(1).meg);
        MEGst = zeros(win, nstim);  % single trial responses
        MEGcm = zeros(win,ncomp);   % ncomp MEGs averaged over nstim
        for icomp = 1:ncomp
            %display( [is icomp])
            MEGcomp = Ysort(is).megcomp(icomp).meg;
            MEGst = MEGst + MEGcompW(icomp)*MEGcomp;
            MEGcm(:,icomp) = MEGcompW(icomp)*mean(MEGcomp,2);
        end
        MEG(is).megst = MEGst; % single trial
        MEG(is).megcm = MEGcm; % components
        MEG(is).megm = mean(MEGst,2); % total mean
    end
end