% maxforce.m - Maximum feasible hand force in each exertion direction.
%
% Reviewers asked what percentage of the model's maximum force the 20, 50 and
% 100 N loads represent. This script answers that exactly, without a search.
%
% Why a single linear program is enough
% -------------------------------------
% Every task uses the same pose, so the muscle moment matrix and the glenoid
% stability matrix are identical across tasks; only the right-hand sides change,
% and both change AFFINELY with the hand force magnitude F:
%
%     beq(F)   = beq_0   + F*dbeq            (required joint moments)
%     bineq(F) = bineq_0 - F*(u*dghjrf)      (stability, load contribution)
%
% Activations and F therefore satisfy one linear system, so the feasible set in
% (a, F) is a polyhedron and its projection onto F is an interval - there can be
% no gaps. Treating F as a decision variable and maximising it gives the exact
% limit in one linprog call per direction and model.
%
% The constraint construction below mirrors feasiblesolns.m exactly: the 1e-5
% moment zeroing, usemuscles = 2:end (dropping the conoid ligament), the reserve
% diagonal, the stability matrix u, the [Fx Fz Fy] column order and the (x,z,y)
% permutation of the load vector. Keep the two in step if either is edited.
%
% Direction naming follows the manuscript, i.e. the direction of the EXERTION,
% which is opposite to the external hand force stored in inputs/tasks.xlsx.
%
% Outputs: printed tables and ../statistics/maxforce.csv

clear all
close all
clc

scriptdir = fileparts(mfilename('fullpath'));
root = fullfile(scriptdir, '..');
outroot = fullfile(root, 'outputs');
statsdir = fullfile(root, 'statistics');
if ~exist(statsdir, 'dir'); mkdir(statsdir); end

% Exertion direction -> the tasks that load it, in ascending magnitude
directions = {'upward', 'downward', 'rightward', 'leftward'};
dirtasks = {{'task05','task06','task07'}, {'task02','task03','task04'}, ...
            {'task11','task12','task13'}, {'task08','task09','task10'}};
magnitudes = [20 50 100];
unloaded = 'task01';

% Model definitions - identical to feasiblesolns.m
select_dofs = {[13:15], [13:17], [10:17], [7:17]};
select_dofs_names = {'3dof', '5dof', '8dof', '11dof'};
reserves_all = [1 1 1 1 1 1, 10 10 10, 10 10 10, 1 1 1, 1 1, 1 1];

u = [0.561, 0, 1;
    0.437, -1, 1;
    0.320, -1, 0;
    0.482, -1, -1;
    0.598, 0, -1;
    0.504, 1, -1;
    0.366, 1, 0;
    0.443, 1, 1];

readtask = @(t, w) readmatrix(fullfile(outroot, t, strcat(t, '-', w, '.csv')));

% Pose-dependent quantities: identical for every task, so read once
musclemoments = readtask(unloaded, 'maxmoments');
musclemoments(abs(musclemoments) < 1e-5) = 0;
[~, nmuscles] = size(musclemoments);
usemuscles = 2:nmuscles;                       % skip conoid ligament
Fiso = readtask(unloaded, 'Fiso');
Fiso = Fiso(usemuscles, usemuscles);
ghjrf_muscle = readtask(unloaded, 'GHJRFmuscle');
ghjrf_muscle = ghjrf_muscle(:, usemuscles);

% Unloaded (gravity-only) right-hand sides
id0 = readtask(unloaded, 'idmoments');
gh0 = readtask(unloaded, 'GHJRFloads');

% Per-newton slopes, least squares through the three magnitudes
dbeq = zeros(length(id0), numel(directions));
dgh = zeros(3, numel(directions));
for d = 1:numel(directions)
    num_id = zeros(size(id0)); num_gh = zeros(size(gh0)); den = 0;
    for k = 1:numel(magnitudes)
        m = magnitudes(k);
        num_id = num_id + m * (readtask(dirtasks{d}{k}, 'idmoments') - id0);
        num_gh = num_gh + m * (readtask(dirtasks{d}{k}, 'GHJRFloads') - gh0);
        den = den + m^2;
    end
    dbeq(:, d) = num_id / den;
    dgh(:, d) = num_gh / den;
end

opts = optimoptions('linprog', 'Display', 'off');

fmt = @(v) regexprep(sprintf('%9.0f', v), 'NaN', 'infeas');

%% Table 1: maximum feasible hand force
fprintf('\nMAXIMUM FEASIBLE HAND FORCE (N), with reserve actuators\n\n');
fprintf('%-11s', 'exertion'); fprintf('%9s', select_dofs_names{:}); fprintf('\n');
Fmax = nan(numel(directions), numel(select_dofs));
for d = 1:numel(directions)
    fprintf('%-11s', directions{d});
    for m = 1:numel(select_dofs)
        Fmax(d, m) = solve_maxforce(select_dofs{m}, musclemoments, usemuscles, Fiso, ...
            ghjrf_muscle, reserves_all, u, id0, gh0, dbeq(:, d), dgh(:, d), 1, true, opts);
        fprintf('%s', fmt(Fmax(d, m)));
    end
    fprintf('\n');
end

%% Table 2: tested loads as a percentage of the directional maximum
fprintf('\n\nTESTED LOADS AS %% OF DIRECTIONAL MAXIMUM\n\n');
fprintf('%-11s%7s%9s', 'exertion', 'model', 'max (N)');
fprintf('%8d N', magnitudes); fprintf('\n');
for d = 1:numel(directions)
    for m = 1:numel(select_dofs)
        fprintf('%-11s%7s%9.0f', directions{d}, select_dofs_names{m}, Fmax(d, m));
        fprintf('%9.0f%%', 100 * magnitudes / Fmax(d, m));
        fprintf('\n');
    end
end

%% Table 3: how much of the limit comes from the reserve actuators
fprintf('\n\nSENSITIVITY TO THE RESERVE ACTUATORS\n');
fprintf('(the 3 and 5 DOF limits are muscle-determined; 8 and 11 DOF are not)\n\n');
fprintf('%-11s%7s%11s%11s%11s%14s\n', 'exertion', 'model', '0.5x res', '1x (used)', '2x res', 'muscles only');
sens = nan(numel(directions), numel(select_dofs), 4);
for d = 1:numel(directions)
    for m = 1:numel(select_dofs)
        v = nan(1, 4);
        scales = [0.5 1.0 2.0];
        for s = 1:3
            v(s) = solve_maxforce(select_dofs{m}, musclemoments, usemuscles, Fiso, ...
                ghjrf_muscle, reserves_all, u, id0, gh0, dbeq(:, d), dgh(:, d), scales(s), true, opts);
        end
        v(4) = solve_maxforce(select_dofs{m}, musclemoments, usemuscles, Fiso, ...
            ghjrf_muscle, reserves_all, u, id0, gh0, dbeq(:, d), dgh(:, d), 1, false, opts);
        sens(d, m, :) = v;
        fprintf('%-11s%7s', directions{d}, select_dofs_names{m});
        fprintf('%11s%11s%11s%14s\n', regexprep(sprintf('%11.0f', v(1)), 'NaN', 'infeas'), ...
            regexprep(sprintf('%11.0f', v(2)), 'NaN', 'infeas'), ...
            regexprep(sprintf('%11.0f', v(3)), 'NaN', 'infeas'), ...
            regexprep(sprintf('%14.0f', v(4)), 'NaN', 'infeasible'));
    end
end

%% Cross-check: does the LP limit agree with which tasks actually solved?
fprintf('\n\nCROSS-CHECK AGAINST THE SAMPLED RESULTS\n');
fprintf('(a task should have a -mineffort- file exactly when its load <= the LP maximum)\n\n');
nagree = 0; ndisagree = 0;
for d = 1:numel(directions)
    for k = 1:numel(magnitudes)
        t = dirtasks{d}{k};
        for m = 1:numel(select_dofs)
            predicted = magnitudes(k) <= Fmax(d, m) + 1e-9;
            solved = isfile(fullfile(outroot, t, 'solutions', ...
                strcat(t, '-mineffort-', select_dofs_names{m}, '.txt')));
            if predicted == solved
                nagree = nagree + 1;
            else
                ndisagree = ndisagree + 1;
                fprintf('  MISMATCH %s %s: LP says %d, files say %d\n', ...
                    t, select_dofs_names{m}, predicted, solved);
            end
        end
    end
end
% the unloaded task should be feasible in every model
for m = 1:numel(select_dofs)
    solved = isfile(fullfile(outroot, unloaded, 'solutions', ...
        strcat(unloaded, '-mineffort-', select_dofs_names{m}, '.txt')));
    if solved; nagree = nagree + 1; else
        ndisagree = ndisagree + 1;
        fprintf('  MISMATCH %s %s: unloaded task has no solution file\n', unloaded, select_dofs_names{m});
    end
end
fprintf('\n  %d of %d task x model combinations agree\n', nagree, nagree + ndisagree);

%% Write the results
rows = {};
for d = 1:numel(directions)
    for m = 1:numel(select_dofs)
        rows(end+1, :) = {directions{d}, select_dofs_names{m}, Fmax(d, m), ...
            100*magnitudes(1)/Fmax(d, m), 100*magnitudes(2)/Fmax(d, m), 100*magnitudes(3)/Fmax(d, m), ...
            sens(d, m, 1), sens(d, m, 3), sens(d, m, 4)};
    end
end
T = cell2table(rows, 'VariableNames', {'Exertion', 'Model', 'MaxForce_N', ...
    'Pct_20N', 'Pct_50N', 'Pct_100N', 'MaxForce_halfReserves', 'MaxForce_doubleReserves', ...
    'MaxForce_musclesOnly'});
writetable(T, fullfile(statsdir, 'maxforce.csv'));
fprintf('\nWritten: %s\n', fullfile(statsdir, 'maxforce.csv'));


% Solve max F for one direction and model.
% resscale: multiplier on the reserve strengths.  usereserves: include them at all.
function [Fmax, exitflag] = solve_maxforce(usedofs, musclemoments, usemuscles, ...
        Fiso, ghjrf_muscle, reserves_all, u, id0, gh0, dbeq_d, dgh_d, resscale, usereserves, opts)

    musclemoments_use = musclemoments(usedofs, usemuscles);
    if usereserves
        res = diag(reserves_all(usedofs) * resscale);
        maxmoments = horzcat(musclemoments_use, res, -res);
    else
        maxmoments = musclemoments_use;
    end
    [ndofs_use, ntorques_use] = size(maxmoments);
    nmuscles_use = ntorques_use - 2 * ndofs_use * usereserves;

    Fx = ghjrf_muscle(1, :)' .* diag(Fiso);
    Fy = ghjrf_muscle(2, :)' .* diag(Fiso);
    Fz = ghjrf_muscle(3, :)' .* diag(Fiso);
    Fx(nmuscles_use + 1:ntorques_use) = 0;
    Fy(nmuscles_use + 1:ntorques_use) = 0;
    Fz(nmuscles_use + 1:ntorques_use) = 0;
    R = [Fx, Fz, Fy];

    Aeq = maxmoments;
    beq = id0(usedofs);
    Aineq = u * R';
    bineq = -u * [gh0(1); gh0(3); gh0(2)];

    n = size(Aeq, 2);
    % F becomes the last decision variable
    Aeq_aug = [Aeq, -dbeq_d(usedofs)];
    Aineq_aug = [Aineq, u * [dgh_d(1); dgh_d(3); dgh_d(2)]];
    f = [zeros(n, 1); -1];                       % maximise F
    lb = [zeros(n, 1); 0];
    ub = [ones(n, 1); Inf];

    [x, ~, exitflag] = linprog(f, Aineq_aug, bineq, Aeq_aug, beq, lb, ub, opts);
    if exitflag == 1
        Fmax = x(end);
    elseif exitflag == -3
        Fmax = Inf;                              % unbounded
    else
        Fmax = NaN;                              % infeasible even unloaded
    end
end
