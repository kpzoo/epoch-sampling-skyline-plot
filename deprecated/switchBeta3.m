% Profile likelihood and iterative estimates for varying beta
clearvars; tic;
clc; close all;

% Assumptions and notes
% - account for non-zero sampling intensities between epochs
% - infer a beta skyline via MLEs
% - have knowledge of when sampling started and ended
% - generalised from betaProfLikTrunc3

% Set shadedplot package
addpath(genpath('/Users/kp10/Documents/MATLAB'));

% Set figure defaults
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);
grey1 = 0.8*ones(1, 3); grey2 = 0.5*ones(1, 3);

% Home and folder to save/load
savetrue = 0;
folload = 'test';
folsave = 'nztest';

%% Simulate sampled-coalescent trees

% Possible trajectories to select
type = 3;
trajNames = {'logis', 'exp', 'steep', 'unif_low', 'unif_high', 'boom', 'cyc', 'bottle', 'mesa'};
trajChoice = trajNames{type};
% Set group size
k = 30;
disp(['Group size, k = ' num2str(k)]);

% Set data source
dataStr = folload;

% Read data generated from phylodyn package in R
thisDir = cd;
cd(['using phylodyn/' dataStr '/' trajChoice '_test']);

% Coalescent and sample times
tcoal = csvread('coaltimes.csv');
tsamp = csvread('samptimes.csv');
% Lineages driving each coalescent
coalLin = csvread('coalLin.csv');
% Trajectory
Nt = csvread('trajy.csv');
t = csvread('trajt.csv');
% Heterochronous tree
tree = phytreeread('tree.txt');

% Samples introduced at each sample time
sampIntro = csvread('sampIntro.csv');
% Sampling constants and split
betaTrue = csvread('beta.csv');
lenSplit = csvread('lensplit.csv');
% Start and end times of sampling epochs
sampStart = csvread('sampEnd1.csv');
sampEnd = csvread('sampEnd2.csv');

cd(thisDir);

% Define num of samples and coalescents
nc = length(tcoal); ns = length(tsamp);

% Combine coalescent and sample times
tLin = sort([tcoal' tsamp']);
len = length(tLin); 

% Get whether a coalescent or sample time (set complements)
isamp = ismember(tLin, tsamp);
icoal = ismember(tLin, tcoal);


%% Heterochronous LTT construction

% Construct LTT, must start with sample
nLin = zeros(size(tLin));
nLin(1) = sampIntro(1);
% Counters for samples and coalescents
c_coal = 0; c_samp = 1;
for j = 2:len
    if isamp(j)
        % Sample event has occurred
        c_samp = c_samp + 1;
        nLin(j) = nLin(j-1) + sampIntro(c_samp);
    else
        % Coalescent event has occurred
        c_coal = c_coal + 1;
        nLin(j) = nLin(j-1) - 1;
    end
end
% Lineages that drive stated events (nLin is after events originally)
nLinPre = nLin;
nLinPre(icoal) = nLinPre(icoal) + 1;
nLinPre(isamp) = nLinPre(isamp) - sampIntro';

% Check have lineage pre-post relationship right
if ~ all(nLinPre(2:end) == nLin(1:end-1))
    error('The nLin and nLinPre relationship is wrong');
end
% Check that all events used and that lineages match
if c_samp ~= ns || c_coal ~= nc || ~all(coalLin' == nLinPre(icoal))
    error('Computed LTT incorrectly');
end

% Times when first get 2 samples and used all samples
id0 = find(nLin == 2, 1, 'first');
idend = find(cumsum(isamp) == ns, 1, 'first');
% Truncat tLin and nLin
t0 = tLin(id0); tend = tLin(idend);
tLinOrig = tLin; nLinOrig = nLin; nLinPreOrig = nLinPre;
tLin = tLin(id0:idend);
nLin = nLin(id0:idend); nLinPre = nLinPre(id0:idend);

% Limit time range of t
idlim = find(t >= t0 & t <= tend);
t = t(idlim); Nt = Nt(idlim);

% Adjust sample epoch counts for start (1 sampling event gone)
sampEpoch = lenSplit;
sampEpoch(1) = sampEpoch(1) - 1;

% Truncate the indicator variables
icoal = icoal(id0:idend); isamp = isamp(id0:idend);
% Redefine ns as less sample events now
ns = length(find(isamp));

% Number of intervals and lengths (times)
dtLin = diff(tLin); numInt = length(dtLin);
% The lineage count over dtLin(i) is nLinPre(i)
lendt = length(dtLin); cumInt = cumsum(dtLin);

% Compute components for grouped skyline on each interval
alpha = 0.5*nLin(1:lendt).*(nLin(1:lendt) - 1);
compCoal = alpha.*dtLin;

% Ids in LTT for various sampling epochs
nepoch = length(lenSplit);
epochtLinID = cell(1, nepoch); epochdel = epochtLinID;
for i = 1:nepoch
    % Time indices in epoch - since starting from last sample in previous
    % epoch use tLin > vs >=
   epochtLinID{i} = find(tLin > sampStart(i) & tLin <= sampEnd(i)); 
   % Interval indices in epoch
   epochdel{i} = find(cumInt > sampStart(i)-t0 & cumInt <= sampEnd(i)-t0);
end
% Remaining interval ids have no sampling intensity (beta = 0 here)
nonzeroID = cell2mat(epochdel);
remdel = setdiff(1:lendt, nonzeroID);
if isempty(remdel)
    disp('No zero sampling intensity intervals');
end

% Test num samples in each epoch
sampEpoch2 = zeros(nepoch, 1);
for i = 1:nepoch
    sampEpoch2(i) = sum(isamp(epochtLinID{i}));
end
if ~all(sampEpoch == sampEpoch2)
    error('Inconsistent sample event sums');
else
    clear('sampEpoch2');
end

% LTT with indicators of sampling epochs, and trace of events
figure;
subplot(3, 1, 1:2); % doubles size of this subplot
stairs(tLin, nLin, 'color', grey1, 'linewidth', 2);
hold on;
for i = 1:nepoch
    ep = epochtLinID{i};
    stairs(tLin(ep), nLin(ep), 'g', 'linewidth', 2);
end
hold off; grid off; box off;
ylabel('LTT');
xlim([t0 tend]);
subplot(3, 1, 3);
stem(tsamp, ones(size(tsamp)), 'c', 'Marker', 'none');
hold on;
stem(tcoal, ones(size(tcoal)), 'color', grey2, 'Marker', 'none');
hold off; grid off; box off;
ylim([0 1.1]);
h = gca; h.YTick = [0 1];
xlim([t0 tend]);
xlabel('time into past');


%% Estimate beta for each epoch with switching of likelihoods

% Number groups and group sizes
krem = rem(lendt, k);
if krem == 0
    % Perfect group divisions
    grpSz = k*ones(1, lendt/k);
else
    % Last group has under k events
    grpSz = [k*ones(1, floor(lendt/k)) krem];
end
nGrp = length(grpSz);

% Indices of groups
idGrp = cell(1, nGrp);
% Length, del, and sum(alpha*del) over each group
A = zeros(1, nGrp); del = A;
% Num of samples/coalescents in each group
s = zeros(1, nGrp); c = s;
% Times at end of each group
tEnd = zeros(1, nGrp);

% Get grouped preferential skyline components
jstart = 1;
for i = 1:nGrp
    % End-indices of group from dtLin
    jstop = jstart + grpSz(i) - 1;
    ids = jstart:jstop;
    idGrp{i} = ids;
    
    % Component sums: A (coalescent) and B (sample)
    A(i) = sum(compCoal(ids));
    % Event types in group
    c(i) = sum(icoal(ids));
    s(i) = sum(isamp(ids));
    if sum([s(i) c(i)]) ~= grpSz(i)
        error('Sum of events inconsistent');
    end
    
    % Group length in time
    del(i) = sum(dtLin(ids));
    % Update end indices
    jstart = jstop + 1;
    
    % End-time of group
    tEnd(i) = tLin(jstop+1);
end

% Skyline under true beta, compSamp checks for beta = 0 regions
[NTrue, compSampTrue] = getPrefSkyBetaEpoch(betaTrue, c, s, nepoch, nGrp,...
    compCoal, dtLin, idGrp, epochdel);
if ~all(compSampTrue(remdel) == 0)
    error('Incorrect account for periods of no sampling');
end

% Num sample events completed
s_ep = zeros(1, nepoch);
for j = 1:nepoch
    s_ep(j) = sum(isamp(epochdel{j}));
end

% Iterative estimate of beta in epochs
[betaEst, NIter, betaIters, dbeta] = mleBetaEpoch(numInt, s, c, nepoch, idGrp,...
    nGrp, compCoal, dtLin, s_ep, epochdel);
betaEst = betaEst'; NIter = NIter'; dbeta = dbeta';
% Max beta derivative
disp(['Worst deriv beta = ' num2str(max(abs(dbeta)))]);


%% Interpolation and plotting

% Vectors for classic grouped skylines
tTrue0 = [t0 tEnd]; tEst0 = [t0 tEnd];
NTrue0 = [NTrue(1) NTrue]; NIter0 = [NIter(1) NIter]; 

% Interpolation to NtZ grid and stats
[NTruet, ~] = getInterp(tTrue0, NTrue0, t, Nt);
[NItert, ~] = getInterp(tEst0, NIter0, t, Nt);

% Fisher information confidence intervals
F = s_ep./(betaEst.^2);
conf = 2./sqrt(F);

% Disaggregrate beta true and estimates to intervals
betaTrue0 = zeros(1, numInt); betaEst0 = betaTrue0;
betaUpper0 = betaTrue0; betaLower0 = betaTrue0;
for j = 1:nepoch
    betaTrue0(epochdel{j}) = betaTrue(j);
    betaEst0(epochdel{j}) = betaEst(j);
    betaUpper0(epochdel{j}) = betaEst(j) + conf(j);
    betaLower0(epochdel{j}) = betaEst(j) - conf(j);
end
betaTrue0 = [betaTrue0(1) betaTrue0];
betaEst0 = [betaEst0(1) betaEst0];
betaUpper0 = [betaUpper0(1) betaUpper0];
betaLower0 = [betaLower0(1) betaLower0];

% Ids with zero sampling rate (epoch ended)
idzEst0 = find(betaEst0 == 0);
idzTrue0 = find(betaTrue0 == 0);
% Ensure zeros match
if ~all(idzEst0 == idzTrue0)
    error('Zero intensity epochs do not match');
else
    idzero = idzEst0;
    clearvars idzEst0 idzTrue0
    idnon = setdiff(1:length(betaEst0), idzero);
end

% Plot all estimated skyline on same figure
figure;
plot(t, Nt, 'k--', 'linewidth', 2);
hold on;
stairs(t, NTruet, 'color', grey2, 'linewidth', 2);
stairs(t, NItert, 'c', 'linewidth', 2);
hold off; grid off; box off;
xlim([t0 tend]);
h1 = gca; h1.YScale = 'log';
ylabel('$\hat{N}$');
xlabel(['time into past, $k = ' num2str(k) '$']);
if savetrue
    cd(folsave);
    saveas(gcf, ['Nep_' num2str(k) '_' num2str(type)], 'fig');
    saveas(gcf, ['Nep_' num2str(k) '_' num2str(type)], 'epsc');
    cd(thisDir);
end

% Plot true and estimated epoch beta
% figure;
% hold on;
% plotErrBnd2(gca, tLin', betaEst0', betaLower0', betaUpper0', 'c', 1);
% stairs(tLin, betaTrue0, '--', 'color', grey2, 'linewidth', 2);
% %stairs(tLin(idzero), zeros(1, length(idzero)), 'ko');
% hold off; grid off; box off;
% xlim([t0 tend]);
% xlabel(['time into past']);
% ylabel('$\beta$');
% if savetrue
%     cd(folsave);
%     saveas(gcf, ['betaep_' num2str(k) '_' num2str(type)], 'fig');
%     saveas(gcf, ['betaep_' num2str(k) '_' num2str(type)], 'epsc');
%     cd(thisDir);
% end

figure;
[ha, hb, hc] = shadedplot(tLin, betaLower0, betaUpper0);
ha(2).FaceAlpha = 1; ha(2).FaceColor = 'c';
hb.Color = [1 1 1]; hc.Color = [1 1 1];
hold on;
stairs(tLin, betaTrue0, '-', 'color', grey1, 'linewidth', 2);
hold off; grid off; box off;
xlim([t0 tend]);
xlabel('time into past');
ylabel('$\beta$');
if savetrue
    cd(folsave);
    saveas(gcf, ['betaep_' num2str(k) '_' num2str(type)], 'fig');
    saveas(gcf, ['betaep_' num2str(k) '_' num2str(type)], 'epsc');
    cd(thisDir);
end

% Combine beta and N skylines 
figure;
subplot(2, 1, 1);
plot(t, Nt, 'k--', 'linewidth', 2);
hold on;
stairs(t, NTruet, 'color', grey2, 'linewidth', 2);
stairs(t, NItert, 'c', 'linewidth', 2);
hold off; grid off; box off;
xlim([t0 tend]);
h1 = gca; h1.YScale = 'log';
ylabel('$\hat{N}$');
xlabel(['$k = ' num2str(k) '$']);
subplot(2, 1, 2);
%hold on;
%plotErrBnd2(gca, tLin', betaEst0', betaLower0', betaUpper0', 'c', 1);
%stairs(tLin, betaTrue0, '--', 'color', grey2, 'linewidth', 2);
[ha, hb, hc] = shadedplot(tLin, betaLower0, betaUpper0);
ha(2).FaceAlpha = 1; ha(2).FaceColor = 'c';
hb.Color = [1 1 1]; hc.Color = [1 1 1];
hold on;
stairs(tLin, betaTrue0, '-', 'color', grey1, 'linewidth', 2);
hold off; grid off; box off;
xlim([t0 tend]); h = gca; h.YLim(1) = 0;
xlabel('time into past');
ylabel('$\hat{\beta}$');
if savetrue
    cd(folsave);
    saveas(gcf, ['skycomb_' num2str(k) '_' num2str(type)], 'fig');
    saveas(gcf, ['skycomb_' num2str(k) '_' num2str(type)], 'epsc');
    cd(thisDir);
end

% Examine convergence of beta
backid = 30;
figure;
plot(betaIters(:, end-backid:end)', 'linewidth', 2);
grid off; box off;
xlabel('$\beta$ iterations');

% Log time
tsim = toc/60;
disp(['Execution time: ' num2str(tsim) ' mins']);