% Computes beta for a given preferential skyline
function [beta, delN] = getBetaEpoch(dtLin, Ngrp, numInt, s_ep, nepoch, idGrp,...
    nGrp, epochdel)

% Assumptions and notes
% - epochs with different betas
% - gives beta conditioned on preferential skyline 
% - used iteratively with getPrefSkyBetaEpoch MLEs 

% Disaggregrate Ngrp to intervals
Nint = zeros(1, numInt);
for i = 1:nGrp
    Nint(idGrp{i}) = Ngrp(i);
end

% Epoch components of beta estimate
beta = zeros(1, nepoch);
delN = zeros(1, nepoch);
for j = 1:nepoch
    % Epoch sum of interval width and population size
    ep = epochdel{j};
    delN(j) = dtLin(ep)*Nint(ep)';
    % MLE of beta from each epoch likelihood
    beta(j) = s_ep(j)/delN(j);
end
   