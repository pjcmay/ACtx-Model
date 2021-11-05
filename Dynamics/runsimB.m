function  [Y, Init] = runsimB(INP, STIMSEQ, Mod, SimCond, Init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run simulation: This version includes adapting afferent synapses
% to one input area, and an explicit Wii for including TRN inhibition of
% MGV
% SimCond.SaveMode: 1 - save only u variable
%                   2 - save all variables
%                   3 - save none (but allows MEGcomp to be harvested) -
%                       NOT IMPLEMENTED YET           
% Copyright Patrick J. C. May January 2016
% (Additional comments: August 2017, Nina Haertwich [N])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y(length(STIMSEQ)).placeholder = [];
options = odeset('RelTol',1e-5,'AbsTol',1e-8,'InitialStep',0.002);
nUnits = Mod.nUnits; % total number of columns [N]

if isempty(Init)
    initU = [zeros(1,nUnits) zeros(1,nUnits) ones(1,nUnits)...
             ones(1,nUnits)]'; % row vector, length 4*(tot # of columns) [N]
    initAaff = ones(Mod.nIn,1);
else
    initU = Init.U;
    initAaff = Init.Aaff;
end

for ii = 1:length(STIMSEQ) % 1:(number of stimuli per SOI) [N]
    
    if SimCond.Display == 1
        display(ii)
        
    end
    
    Faff = INP(STIMSEQ(ii)).inp; % afferent firing rate: inp-matrix made up
    % of row vectors with 'tone modulated by tuning curve embedded in silence'
    % for given column (c.f. makeTonestims.m and addspans.m) [N] 
    
    Faff0 = Faff(1:Mod.nIn,:); % initial condition for afferent firing rate=
    % first column of inp-matrix i.e. Faff-matrix [N]
    
    finIndx = INP(STIMSEQ(ii)).inpspan; % final index = input span, number
    % of entries for silence+tone for single column [N]
    
    Break = INP(STIMSEQ(ii)).break; % baseline [N]
    
    Totspan = (0.001:0.001:finIndx/1000); % speeds up simulations
    
    Points = [0 Break length(Totspan)]; % important points in time: start of
    % 'stimulus-unit'
    
    U = zeros(finIndx, length(initU));
    
    Aaff = zeros(finIndx, length(initAaff));
    
    for ispan = 1:length(Break)+1
        
        p1 = Points(ispan)+1;
        p2 = Points(ispan+1);
        
        if Mod.tau_aff > 0 % if synaptic depression occurs, 
                           % solve differential eqn,
            % note that Aspan is a matrix [N]
            [~, Aspan] = ode45(@(t,y) EqnsAaff(t, y, Mod, Faff0),...
                Totspan(p1:p2), initAaff, options);
            
        else
            Aspan = ones(size(Aaff(p1:p2,:))); % else afferent synaptic
            %strength is constant [N]
        end
        
        % Afferent input
        %[size(Aaff) size(p1:p2) size(Aspan)]
        Aaff(p1:p2,:) = Aspan;
        Syninp = zeros(size(Faff));
        temp = repmat(Aaff',Mod.nInArea,1);
        Syninp(1:Mod.nInArea*Mod.nIn,:) = temp;
        Syninp = Syninp.*Faff;
       
        %Syninp = Faff;
        Inputs.Syninp = Syninp;
        
        % Run equations
        [~, Uspan] = ode45(@(t,y) EqnsBhet(t, y, Mod, Inputs),...
            Totspan(p1:p2), initU, options);
        initU = Uspan(end,:);
        initAaff = Aspan(end,:);
        U(p1:p2,:) = Uspan;
        
    end
    
    % Gather results and save new initial conditions
    if SimCond.DoMEG == 1 && SimCond.SaveMode < 2
        Ytemp.u  = U(:,1:nUnits);
        Ytemp.ui = U(:,nUnits+1:2*nUnits);
        Ytemp.a1 = U(:,2*nUnits+1:3*nUnits);
        Ytemp.a = U(:,3*nUnits+1:4*nUnits); %Note: calcMEGcomp0 accepts Y.a
        Ytemp.aaff = Aaff;
        Ytemp = calcMEGcomp0(Ytemp,Mod,INP,STIMSEQ(ii));
        Y(ii).megcomp = Ytemp.megcomp;
    end        
    switch SimCond.SaveMode
        case 2
        Y(ii).u = U(:,1:nUnits);
        Y(ii).ui = U(:,nUnits+1:2*nUnits);
        Y(ii).a1 = U(:,2*nUnits+1:3*nUnits);
        Y(ii).a = U(:,3*nUnits+1:4*nUnits); %Note: calcMEGcomp0 accepts Y.a
        Y(ii).aaff = Aaff;
        case 1
        Y(ii).u = U(:,1:nUnits);
    end    
end
Init.U = initU;
Init.Aaff = initAaff;
if SimCond.DoMEG == 1 && SimCond.SaveMode == 2
    % Note: this can be run outside this function also
    Y = calcMEGcomp0(Y,Mod,INP,STIMSEQ);
end

