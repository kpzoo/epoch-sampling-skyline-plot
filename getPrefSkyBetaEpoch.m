% Computes preferential skyline for a given beta and k
function [Navg, compSamp] = getPrefSkyBetaEpoch(beta, c, s, nepoch, nGrp,...
    compCoal, dtLin, idGrp, epochdel)

% Assumptions and notes
% - generalised for beta epochs (time-variation)
% - code cleaned, new inputs: s, c, idGrp
% - gives preferential skyline conditioned on beta
% - used iteratively with getBeta MLEs 
% - aim is joint estimation of N(t) and beta

% Input of beta must match epoch number
if length(beta) ~= nepoch
    error('Not enough beta values supplied');
end

% Compute sample components for grouped skyline in each interval 
compSamp = zeros(size(compCoal));
for j = 1:nepoch
    for i = epochdel{j}
        compSamp(i) = beta(j)*dtLin(i);
    end
end

% Main grouped, sampled skyline
Navg = zeros(1, nGrp);

% Construct grouped preferential skyline
for i = 1:nGrp
    % End-indices of group from dtLin
    ids = idGrp{i};
    
    % Component sums: A (coalescent) and B (sample)
    A = sum(compCoal(ids));
    B = sum(compSamp(ids));
    
    % Compute skyline based on coalescents vs samples in group
    if c(i) > s(i)
        % If more coalescent events
        evf = c(i) - s(i);
        l = evf/(2*A);
        % MLE under this condition
        Navg(i) = (l + sqrt(l^2 + B/A))^(-1);
    elseif c(i) < s(i)
        % If more sampling events
        evf = s(i) - c(i);
        l = evf/(2*B);
        % MLE under this condition
        Navg(i) = l + sqrt(l^2 + A/B);
    else
        % Equal numbers of each event type
        Navg(i) = sqrt(A/B);
    end
end

% Check all population sizes valid
if any(isnan(Navg)) || any(isinf(Navg))
    assignin('base','NavgErr', Navg);
    error('Some N values are inadmissible');
end
