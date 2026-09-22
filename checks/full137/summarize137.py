"""Compare the 137-element and 42-element models on the 50 N downward exertion (task03).

Reports the statistics the manuscript uses, so the effect of the muscle-element
reduction can be read directly against the published findings:

  1. mean interquartile activation range by muscle group and DOF model
  2. mean and minimum normalised quadratic effort by DOF model
  3. feasibility, redundancy ordering, and solution-space dimension

The trapezius clavicular attachment is excluded from the range statistics, matching
dataanalysis.R, because it cannot contribute to equilibrium below 11 DOF and its
uniform 0-100% distribution would otherwise dominate the thoracoscapular mean.

Run after inversedynamics137.py and feasiblesolns137.m.  Writes nothing.
"""
import csv, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
TASK = 'task03'
MODELS = {'3dof': 3, '5dof': 5, '8dof': 8, '11dof': 11}
GROUPS = ['thoracoscapular', 'thoracohumeral', 'glenohumeral', 'elbow']

BASE2GROUP = {}
for k in ['trap_scap','trap_clav','lev_scap','pect_min','rhomboid','serr_ant']: BASE2GROUP[k]='thoracoscapular'
for k in ['delt_scap','delt_clav','coracobr','infra','ter_min','ter_maj','supra','subscap']: BASE2GROUP[k]='glenohumeral'
for k in ['lat_dorsi','pect_maj_t','pect_maj_c']: BASE2GROUP[k]='thoracohumeral'

def names_from_model(path):
    return re.findall(r'<Thelen2003Muscle name="([^"]+)"', open(path).read())

def base(n):
    n = re.sub(r'_?\d+$', '', n)                      # full model: trap_scap_3, trap_scap10
    n = re.sub(r'_(S|M|I|P|A|H|U)$', '', n)           # reduced model: trap_scap_S
    return n

def muscle_table(modelfile):
    names = names_from_model(modelfile)[1:]           # drop conoid, matching usemuscles = 2:end
    return names, [BASE2GROUP.get(base(n), 'elbow') for n in names]

ARMS = {
    '42':  (os.path.join(ROOT, 'outputs', TASK, 'solutions'),
            os.path.join(ROOT, 'models', 'shoulder', 'das3_clav_scap_orig_PhD-Thelen-simplified.osim')),
    '137': (os.path.join(HERE, 'outputs', TASK, 'solutions'),
            os.path.join(ROOT, 'models', 'shoulder', 'das3_clav_scap_orig_PhD-Thelen.osim')),
}

def load(arm, model):
    d, mf = ARMS[arm]
    names, groups = muscle_table(mf)
    f = os.path.join(d, f'{TASK}-samples-{model}-all.txt')
    if not os.path.exists(f): return None, names, groups
    return np.loadtxt(f), names, groups

line = lambda c='-', n=76: print(c*n)

print(); line('=')
print(f'137-ELEMENT vs 42-ELEMENT MODEL   |   {TASK}: 50 N downward exertion')
line('=')

# ---- 1. interquartile activation range by muscle group -------------------------
print('\n1. MEAN INTERQUARTILE ACTIVATION RANGE (% maximum), by muscle group')
print('   (trapezius clavicular attachment excluded, as in dataanalysis.R)\n')
print(f"   {'model':7s}{'group':18s}{'42 elements':>14s}{'137 elements':>15s}{'change':>10s}")
iqr = {}
for m in MODELS:
    for g in GROUPS:
        vals = []
        for arm in ('42', '137'):
            S, names, groups = load(arm, m)
            if S is None: vals.append(np.nan); continue
            idx = [i for i, n in enumerate(names)
                   if groups[i] == g and not n.startswith('trap_clav')]
            A = S[:, idx] * 100
            q = np.percentile(A, [25, 75], axis=0)
            vals.append((q[1] - q[0]).mean())
        iqr[(m, g)] = vals
        if np.isnan(vals[0]) and np.isnan(vals[1]): continue
        ch = f'{vals[1]-vals[0]:+.1f}' if np.isfinite(vals[0]) and np.isfinite(vals[1]) else '-'
        f0 = f'{vals[0]:.1f}%' if np.isfinite(vals[0]) else 'infeasible'
        f1 = f'{vals[1]:.1f}%' if np.isfinite(vals[1]) else 'infeasible'
        print(f'   {m:7s}{g:18s}{f0:>14s}{f1:>15s}{ch:>10s}')
    print()

# ---- 2. normalised quadratic effort -------------------------------------------
print('\n2. NORMALISED QUADRATIC EFFORT (sum of squared activations / n actuators)\n')
print(f"   {'model':7s}{'':4s}{'42: mean':>11s}{'137: mean':>12s}{'':4s}{'42: min':>10s}{'137: min':>11s}{'min change':>12s}")
for m in MODELS:
    mean_v, min_v = [], []
    for arm in ('42', '137'):
        d, _ = ARMS[arm]
        c = os.path.join(d, f'{TASK}-cost-{m}-all.txt')
        o = os.path.join(d, f'{TASK}-mineffort-{m}.txt')
        mean_v.append(np.loadtxt(c)[:, 1].mean() if os.path.exists(c) else np.nan)
        if os.path.exists(o):
            x = np.loadtxt(o); min_v.append((x**2).sum()/len(x))
        else: min_v.append(np.nan)
    fm = lambda v: f'{v:.5f}' if np.isfinite(v) else 'infeas'
    ch = f'{100*(min_v[1]/min_v[0]-1):+.1f}%' if np.isfinite(min_v[0]) and np.isfinite(min_v[1]) else '-'
    print(f'   {m:7s}{"":4s}{fm(mean_v[0]):>11s}{fm(mean_v[1]):>12s}{"":4s}{fm(min_v[0]):>10s}{fm(min_v[1]):>11s}{ch:>12s}')

# ---- 3. do the published findings survive? ------------------------------------
print('\n\n3. ARE THE PUBLISHED FINDINGS UNCHANGED?\n')
print(f"   {'model':7s}{'feasible (42 / 137)':>22s}{'redundancy ordering preserved':>32s}")
for m in MODELS:
    fe = []
    for arm in ('42', '137'):
        d, _ = ARMS[arm]
        fe.append(os.path.exists(os.path.join(d, f'{TASK}-mineffort-{m}.txt')))
    ords = []
    for k in (0, 1):
        v = [iqr[(m, g)][k] for g in GROUPS]
        ords.append('n/a' if any(np.isnan(v)) else ('yes' if v == sorted(v) else 'no'))
    # The ordering is only expected once the scapular joints are balanced: below 8 DOF
    # the thoracoscapular muscles meet no constraint and are uniform, so they are the
    # widest rather than the narrowest group. What matters is that both models agree.
    agree = 'both agree: ' + ords[0] if ords[0] == ords[1] else f'DISAGREE 42:{ords[0]} 137:{ords[1]}'
    print(f'   {m:7s}{str(fe[0])+" / "+str(fe[1]):>22s}{agree:>32s}')
print('\n   ordering tested: thoracoscapular < thoracohumeral < glenohumeral < elbow')
print('   expected only for 8 and 11 DOF; below that the thoracoscapular muscles are')
print('   unconstrained and uniformly distributed, so they are the widest group.')

print('\n\n4. SOLUTION-SPACE DIMENSION (muscle activations only)\n')
print(f"   {'model':7s}{'42: n - rank':>15s}{'137: n - rank':>16s}")
for m, nd in MODELS.items():
    print(f'   {m:7s}{42-nd:>15d}{137-nd:>16d}')
print('\n   Each balanced DOF removes exactly one dimension, so the unreduced model')
print('   leaves a far larger space to sample with the same number of samples.\n')
