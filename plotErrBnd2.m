% Plot the prctiles around a mean or median
function plotErrBnd2(currAx, x, y, yL, yU, col, yesstairs)

% Assumptions and notes
% - modified to allow for stairs
% - uses boundedline package and plots on currAx
% - upper and lower bounds yL, yU are absolute (not rel to y)
% - y can be mean, median or any central measure
% - column vectors expected

% Check input structure
if ~iscolumn(x) || ~iscolumn(y) || ~iscolumn(yL) || ~iscolumn(yU)
    error('Expect column vectors - transposing');
end

% If want stairs type plot
if yesstairs
    xorig = x;
    [x, y] = stairs(xorig, y);
    [~, yL] = stairs(xorig, yL);
    [~, yU] = stairs(xorig, yU);
end

% Bounded line with formatting and outline
axes(currAx);
[l, p] = boundedline(x, y, [y - yL yU - y], '-');
p.FaceAlpha = 0.25; p.FaceColor = col;
l.Color = col; l.LineWidth = 2;
outlinebounds(l, p);