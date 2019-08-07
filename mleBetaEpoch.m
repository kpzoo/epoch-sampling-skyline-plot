% Iteratively estimate beta and N simultaneously
function [betaEst, NEst, beta, dbeta] = mleBetaEpoch(numInt, s, c, nepoch, ...
    idGrp, nGrp, compCoal, dtLin, s_ep, epochdel)

% Assumptions and notes
% - have nepoch betas to estimate
% - generalisation of skyMLEBetaFnLik3 for epochs
% - expects truncated t and Nt etc

% Initialise a beta epoch set at random
M = 10^4; beta = zeros(nepoch, M);
beta(:, 1) = unifrnd(10^(-6), 10^2, [nepoch 1]);

% Population groups and epoch del*N
Niter = zeros(nGrp, M-1);
delN = zeros(nepoch, M);

% Iterate until beta converges
for i = 2:M
    % Calculate N(t) mean estimate given beta
    [Niter(1:nGrp, i-1), ~] = getPrefSkyBetaEpoch(beta(:, i-1), c, s, nepoch,... 
        nGrp, compCoal, dtLin, idGrp, epochdel);
    % Epoch beta estimates given Niter
    [beta(:, i), delN(:, i)] = getBetaEpoch(dtLin, Niter(1:nGrp, i-1),...
        numInt, s_ep, nepoch, idGrp, nGrp, epochdel);
end
% Final estimates from iterations
betaEst = beta(:, end);
NEst = Niter(:, end);

% Derivatives of likelihood, should be close to 0
dbeta = zeros(1, nepoch);
for i = 1:nepoch
   dbeta(i) = s_ep(i)/betaEst(i) -  delN(i, end);
end

% % Display results
% disp(['True beta: ' num2str(betaTrue')]);
% disp(['Iter beta: ' num2str(betaEst')]);
% disp(['Derivs across beta epochs: ' num2str(dbeta)]);

