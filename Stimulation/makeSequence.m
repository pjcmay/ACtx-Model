function Seq0 = makeSequence(n,p,mode)
% Builds oddball stimulus sequence
% n is number of stimuli
% p is the probability of the deviant
% Copyright Patrick May, Lancaster University, 2021

if nargin == 2
    mode = 1;
end
if length(p) == 1
    if mode == 1
        % This version ensures that deviants are always separated 
        % by at least one presentation of the standard. Slightly clumsily,
        % the length of the sequence is a multiple of len = round(1/p) 
        len = round(1/p);
        m = round(n/len);
        perm = randperm(len);
        Seq = ones(1,len);
        Seq(perm == len) = 2;
        for i = 1:m-1
            perm = randperm(len);
            block = ones(1,len);
            block(perm == len) = 2;
            if block(1) + Seq(end) == 4  % make sure dev dev does not occur
                perm = randperm(len-1);
                block = ones(1,len-1);
                block(perm == len-1) = 2;
                block = [1 block]; %#ok<AGROW>
            end
            Seq = [Seq block]; %#ok<AGROW>
        end
        Seq0 = Seq;
    elseif mode == 2
        % No dev-dev restrictions
        temp = ones(1,n);
        temp2 = randperm(n);
        temp(temp2(1:ceil(p*n))) = 2;
        Seq0 = temp;
    end
else
    temp = ones(1,n);
    temp2 = randperm(n);
    p0 = 1;
    for i = 1:length(p)
        temp(temp2(p0:p0+(p(i)*n)-1)) = i+1;
        p0 = p0+p*n;
    end
    Seq0 = temp;
end