"""Inverse dynamics for the UNREDUCED 137-element model, all 13 tasks.

Identical to ../../inversedynamics.py except for the model file and the output
location. Requires OpenSim 4.4: 4.5.2 rejects this model at initSystem() because
conoid_lig has a single-point path and 4.5.2 enforces at least two PathPoints.
Inverse dynamics does not involve muscles, so the joint moments it produces are
bit-identical to the 42-element model; only the muscle moment arms, maximum
isometric forces and GH force directions differ.
"""

import sys
from pathlib import Path

# utils.py lives at the project root, two levels up from checks/full137/
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from utils import *


# Define working directory
# Resolve the project root from this script's location (checks/full137/ -> root)
wdir = Path(__file__).resolve().parents[2]
outroot = Path(__file__).resolve().parent / 'outputs'

# Model files
modelname = 'das3_clav_scap_orig_PhD-Thelen.osim'   # full 137-element model
modelname_nomuscle = 'das3_clav_scap_orig_PhD-Thelen-nomuscles.osim'
modelfile = str(wdir / ('models/shoulder/' + modelname))
modelfile_nomuscle = str(wdir / ('models/shoulder/' + modelname_nomuscle))

# Setup files 
genericmotfile_pose = str(wdir / 'inputs/genericpose.mot')
genericmotfile_loads = str(wdir / 'inputs/genericloads.mot')
id_setupfile = str(wdir / 'inputs/idsetup.xml')
so_setupfile = str(wdir / 'inputs/sosetup.xml')
jra_setupfile = str(wdir / 'inputs/jrasetup.xml')
reserves_setupfile = str(wdir / 'inputs/unlimitedreserves.xml')
    
# Task inputs
taskfile = str(wdir / 'inputs/tasks.xlsx')
taskfile = pd.read_excel(taskfile)
# One representative task: 50 N downward exertion
TASK = 'task03'
taskfile = taskfile[taskfile['task'] == TASK].reset_index(drop=True)
ntasks = len(taskfile)

# Run through each task
for task_in in range(ntasks):

    # Get task inputs
    # Force directions: X+ = right; Y+ = up; Z+ = backwards (these are external forces)
    task = taskfile["task"].loc[task_in]
    pose = taskfile.loc[:, 'TH_x':'wrist_flexion'].iloc[task_in].values
    forces = taskfile.loc[:, 'Fx':'Fz'].iloc[task_in].values
    torques = taskfile.loc[:, 'Tx':'Tz'].iloc[task_in].values

    # Check that trunk isn't rotated, otherwise forces need to be rotated.
    if np.any(pose[0:3] != 0) and np.any(forces != 0):
        raise Exception("Trunk is rotated, forces need to be rotated.")

    # Output folder and files
    outfolder = str(outroot / task)
    os.makedirs(outfolder, exist_ok=True)
    outmotfile_pose = outfolder + '/' + task + '-pose.mot'
    outmotfile_loads = outfolder + '/' + task + '-loads.mot'

    # Load model, initialize, and identify relevant model parameters
    model = osim.Model(modelfile)
    state = model.initSystem()
    nDofs = model.getCoordinateSet().getSize()
    handcom = np.zeros(3)
    for axis in range(3):
        handcom[axis] = model.get_BodySet().get('hand_r').get_mass_center().get(axis)

    # Create .mot files for pose and forces
    createmot_staticpose(pose, genericmotfile_pose, outmotfile_pose)
    createmotfile_forces(forces, handcom, torques, genericmotfile_loads, outmotfile_loads)

    # Pose model (need to first convert rotational dofs to radians).
    # Work on a copy: pose is written to the .mot file in degrees above, so mutating
    # it in place would silently corrupt that file if the order were ever changed.
    translationdofs = [3, 4, 5]
    rotationaldofs = np.full(len(pose), True, dtype=bool)
    rotationaldofs[translationdofs] = False
    pose_rad = np.array(pose, dtype=float)
    pose_rad[rotationaldofs] = np.deg2rad(pose_rad[rotationaldofs])
    for dof in range(nDofs):
        coordinate = model.getCoordinateSet().get(dof)
        coordinate.setValue(state, pose_rad[dof])

    # Obtain muscle information (r = moment arms, Fiso = max isometric force, maxmoment = max isometric moment)
    # This will include the conoid ligament.
    # These depend only on the model and the pose, not on the applied load, so they are
    # computed once and reused for every task that shares the first task's pose.
    if task_in == 0:
        Fiso, maxmoment = osim_muscleinfo(model, state)
        pose_reference = np.array(pose, dtype=float)
    elif np.array_equal(np.array(pose, dtype=float), pose_reference):
        pass  # same pose, reuse Fiso and maxmoment from the first task
    else:
        Fiso, maxmoment = osim_muscleinfo(model, state)

    # Obtain input kinematics
    motData = osim.Storage(outmotfile_pose)
    initial_time = motData.getFirstTime()
    final_time = motData.getLastTime()

    # Setting the hand loads
    loads_setupfile = str(wdir / 'inputs/genericloads_setup.xml')
    LoadsTool = osim.ExternalLoads(loads_setupfile, True)
    LoadsTool.setDataFileName(outmotfile_loads)
    LoadsTool.get('externalforce').set_data_source_name(os.path.basename(outmotfile_loads))
    loadssetupfile_task = os.path.basename(outmotfile_pose)[:-9] + '-loadssetup' + '.xml'
    LoadsTool.printToXML(loadssetupfile_task)
    shutil.move(os.path.abspath(loadssetupfile_task), outfolder + '/' + loadssetupfile_task)

    # Run inverse dynamics
    osim_run_id(model=model, id_setupfile=id_setupfile,
        posefile=outmotfile_pose, loadsfile=loadssetupfile_task, ti=initial_time, tf=final_time, resultsdir=outfolder)

    # Read the .sto file output and grab moments (static pose, so only use first row and skip time column)
    moments = readmotfile(outfolder + '/' + os.path.basename(outmotfile_pose)[:-8] + 'id.sto', skiprows=7,
                          dataarray=True)
    moments = moments[0, 1:]
    dofnames = ['TH_x', 'TH_y', 'TH_z', 'TH_transx', 'TH_transy', 'TH_transz',
                'SC_y', 'SC_z', 'SC_x', 'AC_y', 'AC_z', 'AC_x',
                'GH_y', 'GH_z', 'GH_yy', 'EL_x', 'PS_y',
                'wrist_deviation', 'wrist_flexion']
    formatted_moments = ', '.join([f'{dof}: {moment:.3f}' for dof, moment in zip(dofnames, moments)])
    print(formatted_moments)

    # Define scapular coordinate system
    R_scap, Rinv_scap, origin = scap_rotation(model, state)

    # Calculate muscle vector components for GH JRF
    gh_jrf_muscle, allmusclenames = calc_ghjrf_muscle(model=model, state=state, Rinv_scap=Rinv_scap, origin=origin)

    # Perform JRA to get GH joint reaction forces
    osim_run_jra(modelfile=modelfile_nomuscle, so_setupfile=so_setupfile, jra_setupfile=jra_setupfile,
        posefile=outmotfile_pose, loadsfile=loadssetupfile_task, ti=initial_time, tf=final_time, resultsdir=outfolder,
        forcesetfile=reserves_setupfile)

    # GH JRF loads (gravity + external forces) vector components (these are external reaction forces)
    gh_jrf_load = readmotfile(outmotfile_pose[:-9] + '_JointReaction_ReactionLoads.sto', skiprows=12, dataarray=True)
    gh_jrf_load = -gh_jrf_load[0, 64:67]  # OpenSim gives internal, multiplying to get external
    gh_jrf_load = Rx.dot(Rinv_scap.dot(gh_jrf_load))  # Transforming to scap coordinate system

    # Data visualization
    visualfolder = outfolder + '/figs-modelparameters'
    visualpath = visualfolder + '/' + task
    os.makedirs(visualfolder, exist_ok=True)
    visualizeoutputs(moments=moments, dofnames=dofnames, maxforces=Fiso, musclenames=allmusclenames, 
        maxmoments=maxmoment, gh_jrfs=gh_jrf_muscle, savepath=visualpath)

    # Save outputs
    np.savetxt(outfolder + '/' + task + '-maxmoments.csv', maxmoment, delimiter=',') # Max muscle moments
    np.savetxt(outfolder + '/' + task + '-idmoments.csv', moments, delimiter=',') # ID moments
    np.savetxt(outfolder + '/' + task + '-Fiso.csv', Fiso, delimiter=',') # Max muscle forces
    np.savetxt(outfolder + '/' + task + '-GHJRFmuscle.csv', gh_jrf_muscle, delimiter=',') # GH JRF muscle vector components
    np.savetxt(outfolder + '/' + task + '-GHJRFloads.csv', gh_jrf_load, delimiter=',') # GH JRF loads vector components
