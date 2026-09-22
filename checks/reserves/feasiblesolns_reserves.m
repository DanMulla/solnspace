% feasiblesolns_reserves.m - reserve actuator sensitivity check.
%
% Repeats the feasible solution sampling for the 50 N downward exertion (task03)
% under three reserve actuator settings, reproducing the sampling behind Fig D in
% S1 File and the sensitivity statements in the Discussion:
%
%   SCAC5   sternoclavicular and acromioclavicular reserves reduced 10 -> 5 Nm
%   SCAC20  sternoclavicular and acromioclavicular reserves raised  10 -> 20 Nm
%   GH01    glenohumeral reserves reduced 1 -> 0.1 Nm (SC and AC left at 10 Nm)
%
% Everything else (constraints, stability matrix, sample counts, seeds) is
% unchanged from ../../feasiblesolns.m, and the seed counter is reset for each
% variant so differences are attributable to reserve strength alone.
%
% Runtime: roughly 6 min for three variants x four DOF models.

clear all
close all
clc

tic

% Number of points to sample
nsamples = 10^4;

% Set seed for reproducibility
rngseed = 2024;
seedcounter = 0;

% Add to path: the sampler lives at the project root, two levels up
scriptdir = fileparts(mfilename('fullpath'));
root = fullfile(scriptdir, '..', '..');
addpath(genpath(fullfile(root, 'PolytopeSamplerMatlab')))
addpath(genpath(scriptdir))

% Reserve actuator strength does not affect inverse dynamics, so model inputs are
% read from the main results folder; only sampled solutions are written here.
task = 'task03';
outfolder = fullfile(root, 'outputs', task);

% Input files
musclemoments = readmatrix(fullfile(outfolder, strcat(task, '-maxmoments.csv')));
jointmoments = readmatrix(fullfile(outfolder, strcat(task, '-idmoments.csv')));
[ndofs, nmuscles] = size(musclemoments);

% Set tiny muscle moments as 0
musclemoments(abs(musclemoments) < 1e-5) = 0;

% Reserve actuator variants. Order of DOF: trunk (6), SC (3), AC (3), GH (3),
% ELx-PSy (2), wrist (2). The main analysis uses 10 Nm at SC and AC, 1 Nm elsewhere.
variants = {'SCAC5', 'SCAC20', 'GH01'};
variant_strengths = { ...
    [1 1 1 1 1 1,  5  5  5,  5  5  5,  1   1   1,   1 1, 1 1], ...   % SC and AC reduced to 5 Nm
    [1 1 1 1 1 1, 20 20 20, 20 20 20,  1   1   1,   1 1, 1 1], ...   % SC and AC raised to 20 Nm
    [1 1 1 1 1 1, 10 10 10, 10 10 10, 0.1 0.1 0.1,  1 1, 1 1]};      % GH reduced to 0.1 Nm
usereserves = true;

% Select DOFs to balance (Gh, Gh+Elb, Gh+Elb+Ac, Gh+Elb+Ac+Sc)
select_dofs = {[13:15], [13:17], [10:17], [7:17]};
select_dofs_names = {'3dof', '5dof', '8dof', '11dof'};

% Which muscles to use
usemuscles = [2:nmuscles];  % Skip conoid ligament

% GHJRF muscle vectors and external loads
Fiso = readmatrix(fullfile(outfolder, strcat(task, '-Fiso.csv')));
ghjrf_load = readmatrix(fullfile(outfolder, strcat(task, '-GHJRFloads.csv')));
ghjrf_muscle = readmatrix(fullfile(outfolder, strcat(task, '-GHJRFmuscle.csv')));
Fiso = Fiso(usemuscles, usemuscles);
ghjrf_muscle = ghjrf_muscle(:, usemuscles);

% Indicator of issues
unfeasibleprobs = 0;
maxeffort_warnings = 0;

%% Running through each reserve actuator variant
for v = 1:numel(variants)

reserves = diag(variant_strengths{v});
resultsdir = fullfile(scriptdir, 'outputs', variants{v}, task, 'solutions');
if ~exist(resultsdir, 'dir'); mkdir(resultsdir); end
seedcounter = 0;   % reset so every variant draws the same seeds
disp(['=== reserve variant: ' variants{v} ' ===']);

%% Running through different versions of the model
for model_dof = 1:length(select_dofs)

    disp(strcat('Running simulations for:', select_dofs_names{model_dof}));

    %% Subset the data for select dofs and muscles
	usedofs = select_dofs{model_dof};
	musclemoments_use = musclemoments(usedofs, usemuscles);
	jointmoments_use = jointmoments(usedofs);
	reserves_use = reserves(usedofs, usedofs);
	maxmoments = horzcat(musclemoments_use, reserves_use, -reserves_use);
	
	[ndofs_use, ntorques_use] = size(maxmoments);
	nmuscles_use = ntorques_use - 2*ndofs_use;

    %% Defining GHJRF constraints
    % The maximum force each muscle can generate across each axis
    % +Fx: Right; +Fy: Up; +Fz: Backwards
    Fx = ghjrf_muscle(1, :)' .* diag(Fiso);
    Fy = ghjrf_muscle(2, :)' .* diag(Fiso);
    Fz = ghjrf_muscle(3, :)'.* diag(Fiso);
    Fx(nmuscles_use + 1: ntorques_use) = 0;  % 0's for reserves
    Fy(nmuscles_use + 1: ntorques_use) = 0;
    Fz(nmuscles_use + 1: ntorques_use) = 0;

    % Dislocation thresholds: u*R <= 0 (from Halder et al. 2001)
    % 8 linear inequality constraints - Dickerson et al. 2007
    % See: Blache et al. (2016). Sports Biomechanics. 16(1):127-142.
    % Starting from top and going clockwise if looking into glenoid.
    % Note: The z-axis is flipped (forwards is negative) 
    % The stability cone applies to the GHJRF:
    % u * (R'*a + f_load) <= 0
    % (u*R')*a <= -u*f_load
    u = [0.561, 0, 1;
        0.437, -1, 1;
        0.320, -1, 0;
        0.482, -1, -1;
        0.598, 0, -1;
        0.504, 1, -1;
        0.366, 1, 0;
        0.443, 1, 1];
    R = [Fx, Fz, Fy];

    %% Linear problem
    Aeq = maxmoments;
    beq = jointmoments_use;
    Aineq = u * R';
    bineq = -u * [ghjrf_load(1); ghjrf_load(3); ghjrf_load(2)];
    n = size(Aeq, 2);
    lb = zeros(n, 1);
    ub = ones(n, 1);

    %% Optimal solutions (mineffort^2)
    % Using fmincon
    % objective = @(x) 0.5 * x' * (2 * eye(n)) * x;
    objective = @(x) sum(x.^2);
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');
    [xmin, fval_min, exitflag_min, output_min] = fmincon(objective, zeros(n, 1), Aineq, bineq, Aeq, beq, lb, ub, [], options);
    [xmax, fval_max, exitflag_max, output_max] = fmincon(@(x) -objective(x), zeros(n, 1), Aineq, bineq, Aeq, beq, lb, ub, [], options);

    % A positive exitflag alone does not guarantee the returned point satisfies the
    % constraints, so check the residuals explicitly. maxviolation returns the worst
    % violation across the equality, inequality and bound constraints.
    maxviolation = @(x) max([max(abs(Aeq*x - beq)); Aineq*x - bineq; lb - x; x - ub]);
    feastol = 1e-6;
    % Feasibility is decided by the constraint residual, not by the exitflag. fmincon
    % returns exitflag 0 when it hits its iteration limit, which says nothing about
    % whether the point it returned is feasible: for the unloaded 11 DOF problem it
    % returns exitflag 0 with a residual of 1e-15 and an effort value unchanged by
    % raising MaxIterations, while for the 50 N upward 11 DOF problem it returns
    % exitflag 0 with a residual of 2.6. The residual separates these two cases; the
    % exitflag does not, and treating it as the test discards valid solutions.
    minfeasible = maxviolation(xmin) < feastol;
    minconverged = exitflag_min > 0;
    maxfeasible = maxviolation(xmax) < feastol;   % usable as an effort interval bound
    maxconverged = exitflag_max > 0;
    skipeffort = false;
    fprintf('Constraint violation: xmin %.3e, xmax %.3e (tolerance %.0e)\n', ...
        maxviolation(xmin), maxviolation(xmax), feastol);

    if minfeasible && maxfeasible && minconverged && maxconverged
        disp('A feasible solution exists.');
    elseif minfeasible && maxfeasible
        % One or both solves stopped on the iteration limit, but the returned points
        % satisfy the constraints, so they remain valid. Proceed and report which.
        if ~minconverged
            disp('Minimum effort solution is feasible but fmincon stopped on its');
            disp('iteration limit; verify the effort value is stable if in doubt.');
        end
        if ~maxconverged
            disp('Maximum effort solution is feasible but did not converge cleanly,');
            disp('so the sampled effort range may be conservative.');
            maxeffort_warnings = maxeffort_warnings + 1;
        end
    elseif minfeasible
        disp('Maximum effort solution is NOT feasible to tolerance; the effort interval');
        disp('bounds would be invalid, so effort landscape sampling is skipped.');
        maxeffort_warnings = maxeffort_warnings + 1;
        skipeffort = true;
    else
        disp('Failed to find minimum and maximum feasible solutions.');
        disp(output_min.message);
        disp(output_max.message);
        unfeasibleprobs = unfeasibleprobs + 1;
        continue
    end

    writematrix(xmin, fullfile(resultsdir, strcat(task, '-mineffort-', select_dofs_names{model_dof})), 'Delimiter', 'tab');

    %% Feasible solutions (entire space)
    P = struct;
    P.Aeq = Aeq;
    P.beq = beq;
    P.Aineq = Aineq;
    P.bineq = bineq;
    P.lb = lb;
    P.ub = ub;
    opts = default_options();
    seedcounter = seedcounter + 1;
    opts.seed = rngseed + seedcounter;
    o = sample(P, nsamples, opts);
    
    % Convergence and verification steps
    o.summary;
    nsamplesobt = length(o.samples);
    disp('Number of samples outputted:');
    disp(nsamplesobt);
    error = P.Aeq*o.samples - jointmoments_use;
    disp('Max aboslute error:')
    maxerror = max(abs(error), [], 2);
    disp(maxerror);
    distribution_test(o, struct('toPlot', true));
    drawnow()
    saveas(gcf, fullfile(resultsdir, strcat(select_dofs_names{model_dof}, '-uniformity-all.tif')))

    %% Feasible solutions (sampling across effort landscape)
    if skipeffort
        disp('Skipping effort landscape sampling (infeasible maximum effort solution).');
        close all
        continue
    end
    mineffort = sum(xmin);
    maxeffort = sum(xmax);
    intervals = 100;
    effort_intervals = linspace(mineffort, maxeffort, intervals+1);
    effortsolns = zeros(ntorques_use, nsamples);  % Placeholder for storing all solutions
    i = 1;
    for interval = 1:intervals
        P.Aineq = [Aineq; ones(1, n); -ones(1, n)];
        P.bineq = [bineq; effort_intervals(interval+1); -effort_intervals(interval)];

        seedcounter = seedcounter + 1;
        opts.seed = rngseed + seedcounter;
        o_effort = sample(P, nsamples / intervals, opts);

        % Convergence and verification steps
        o_effort.summary;
        nsamplesobt = size(o_effort.samples, 2);

        % Concatenate data
        effortsolns(:, i:i+nsamplesobt-1) = o_effort.samples;
        i = i + nsamplesobt;
    end

    % The sampler returns approximately, not exactly, nsamples/intervals per interval,
    % so the placeholder may be under- or over-filled. Columns 1:i-1 are the ones
    % actually written, so trim to those (length() on a matrix returns its largest
    % dimension, which is not the column count and never triggered the old guard).
    effortsolns = effortsolns(:, 1:i-1);
    % Placing back in o_effort to keep rest of code same
    o_effort.samples = effortsolns;
    error = P.Aeq*o_effort.samples - jointmoments_use;
    disp('Max aboslute error:')
    maxerror = max(abs(error), [], 2);
    disp(maxerror);

    o.rngseed = rngseed;
    o_effort.rngseed = rngseed;
    o_effort.mineffort = xmin;
    o_effort.maxeffort = xmax;

    %% Objectives and constraints

    % All samples
    o.costs(:, 1) = effort(o.samples, 1)/ntorques_use;
    o.costs(:, 2) = effort(o.samples, 2)/ntorques_use;
    o.costs(:, 3:8) = ghjrf(Fiso, o.samples, ghjrf_muscle, ghjrf_load);

    % Low effort solutions
    o_effort.costs(:, 1) = effort(o_effort.samples, 1)/ntorques_use;
    o_effort.costs(:, 2) = effort(o_effort.samples, 2)/ntorques_use;
    o_effort.costs(:, 3:8) = ghjrf(Fiso, o_effort.samples, ghjrf_muscle, ghjrf_load);

    %% Data Visualization

    % Activation landscape
    subplot(1, 2, 1)
    boxplot(o.samples')
    ylim([0 1])
    subplot(1, 2, 2)
    boxplot(o_effort.samples')
    ylim([0 1])
    figureHandle = gcf;
    set(figureHandle, 'Position', [100, 100, 1200, 400]);
    saveas(gcf, fullfile(resultsdir, strcat(select_dofs_names{model_dof}, '-activations.tif')))
    
    % Effort landscape
    subplot(2, 3, 1)
    histogram(o.costs(:, 1))
    xline(mineffort/ntorques_use, 'r', 'LineWidth', 2);
    xline(maxeffort/ntorques_use, 'r', 'LineWidth', 2);
    xlim([0 1])
    subplot(2, 3, 4)
    histogram(o.costs(:, 2))
    xline(sum(xmin.^2)/ntorques_use, 'r', 'LineWidth', 2);
    xline(sum(xmax.^2)/ntorques_use, 'r', 'LineWidth', 2);
    xlim([0 1])
    subplot(2, 3, 2)
    histogram(o_effort.costs(:, 1))
    xline(mineffort/ntorques_use, 'r', 'LineWidth', 2);
    xline(maxeffort/ntorques_use, 'r', 'LineWidth', 2);
    xlim([0 1])
    subplot(2, 3, 5)
    histogram(o_effort.costs(:, 2))
    xline(sum(xmin.^2)/ntorques_use, 'r', 'LineWidth', 2);
    xline(sum(xmax.^2)/ntorques_use, 'r', 'LineWidth', 2);
    xlim([0 1])
    subplot(2, 3, 3)
    boxplot(o_effort.costs(:,1))
    hold on
    plot(1, mineffort/ntorques_use, '-ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    plot(1, maxeffort/ntorques_use, '-ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    ylim([0 1])
    subplot(2, 3, 6)
    boxplot(o_effort.costs(:,2))
    hold on
    plot(1, sum(xmin.^2)/ntorques_use, '-ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    plot(1, sum(xmax.^2)/ntorques_use, '-ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    ylim([0 1])
    saveas(gcf, fullfile(resultsdir, strcat(select_dofs_names{model_dof}, '-effort.tif')))

    %% Save Data
    outfile_act = strcat(task, '-samples-', select_dofs_names{model_dof});
    save(fullfile(resultsdir, strcat(outfile_act, '-all.mat')), 'o');
    save(fullfile(resultsdir, strcat(outfile_act, '-loweffort.mat')), 'o_effort');
    writematrix(o.samples', fullfile(resultsdir, strcat(outfile_act, '-all.txt')), 'Delimiter', 'tab');
    writematrix(o_effort.samples', fullfile(resultsdir, strcat(outfile_act, '-loweffort.txt')), 'Delimiter', 'tab');
    outfile_cost = strcat(task, '-cost-', select_dofs_names{model_dof});
    writematrix(o.costs, fullfile(resultsdir, strcat(outfile_cost, '-all.txt')), 'Delimiter', 'tab');
    writematrix(o_effort.costs, fullfile(resultsdir, strcat(outfile_cost, '-loweffort.txt')), 'Delimiter', 'tab');

    close all
end
end  % reserve variant


disp('Warnings about max effort optimal solution:');
disp(maxeffort_warnings);
disp('Number of unfeasible problems:');
disp(unfeasibleprobs);
toc

%% Functions

function J = effort(activation, power)
    J = sum(activation.^power, 1);
end

function G = ghjrf(fmax, activation, ghjrf_musclevector, ghjrf_arm)
    nmusc = size(fmax, 1);
    fm = activation(1:nmusc, :)' * fmax;  % Muscle forces
    ghjrf_vector = fm * ghjrf_musclevector' + ghjrf_arm';  % GH JRF vector components
    ghjrf = sum(ghjrf_vector.^2, 2).^0.5;  % Resultant
    compression = abs(ghjrf_vector(:, 1));  % glenoid-normal component
    sy = ghjrf_vector(:, 2) ./ compression;  % stability ratio up/down
    sz = ghjrf_vector(:, 3) ./ compression;  % Note: + = posterior
    sy(compression == 0) = 0;  % avoid Inf/NaN when there is no compression
    sz(compression == 0) = 0;  % (matches nan_to_num in utils.py visualizeoutputs)
    G = [ghjrf_vector, ghjrf, sy, sz];
end