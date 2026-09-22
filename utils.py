import csv
import matplotlib.pyplot as plt
import numpy as np
import opensim as osim
import os
import pandas as pd
import shutil


# Defining a rotation around x-axis by 19 degrees for glenoid torsion (Matsumura et al. (2018). JESE Open Access 2:194-199)
theta = np.radians(19)
Rx = np.array([[1, 0, 0],
    [0, np.cos(theta), -np.sin(theta)],
    [0, np.sin(theta), np.cos(theta)]])


def calc_ghjrf_muscle(model, state, Rinv_scap, origin):

    """
    Calculate GH JRF vector components for each muscle for given state.

    """
    
    nmuscles = model.getForceSet().getSize()
    gh_jrf_muscle = np.zeros([3, nmuscles])
    allmusclenames = []
    
    for muscle in range(0, nmuscles):

        activemuscle = model.getMuscles().get(muscle)
        musclename = activemuscle.getName()
        allmusclenames.append(musclename)
        musclepath = activemuscle.getGeometryPath().getCurrentPath(state)
        nattachments = musclepath.getSize()
        muscleattachments = np.zeros([nattachments, 3])
        musclebodies = []

        # Grab muscle attachmenst and transform in scapular coordinate system
        for point in range(nattachments):
            attachment_temp = activemuscle.getGeometryPath().getCurrentPath(state).get(point).getLocationInGround(state)
            attachment_ground = np.array([attachment_temp[0], attachment_temp[1], attachment_temp[2]])
            attachment_scap = Rx.dot(Rinv_scap.dot(attachment_ground - origin))
            muscleattachments[point, :] = np.array([attachment_scap[0], attachment_scap[1], attachment_scap[2]])
            musclebodies.append(musclepath.get(point).getBody().getName())

        # Looping through attachments and identifying humeral attachments
        # Getting the net vector for each muscle humeral attachment (if not attaching on humerus, assigned [0, 0, 0])
        if 'humerus_r' in musclebodies:
            hum_ids = [id for id, body in enumerate(musclebodies) if body == 'humerus_r']
            nattachments_hum = len(hum_ids)
            nethumeralvector = np.zeros([nattachments_hum, 3])

            # Getting the vector from the given humeral attachment relative to point before and after
            # Ignoring forearm muscles
            for i, hum_point in enumerate(hum_ids):

                # Point afterwards - current point
                if hum_point < nattachments - 1 and musclebodies[hum_point + 1] not in ['ulna_r', 'radius_r']:
                    vector1 = muscleattachments[hum_point + 1] - muscleattachments[hum_point]
                    vector1 = vector1 / np.linalg.norm(vector1)
                else:
                    vector1 = np.zeros(3)

                # Point previously - current point
                if hum_point > 0 and musclebodies[hum_point - 1] not in ['ulna_r', 'radius_r']:
                    vector2 = muscleattachments[hum_point - 1] - muscleattachments[hum_point]
                    vector2 = vector2 / np.linalg.norm(vector2)
                else:
                    vector2 = np.zeros(3)

                nethumeralvector[i, :] = vector1 + vector2

            # Summing up the humeral vector components
            nethumeralvector = nethumeralvector.sum(axis=0)

        # Treating these muscles separately as they go from scapula directly to forearm, so indirect effect
        # on GH JRF: this is to match OpenSim's JRA tool
        elif musclename in ['bic_b', 'tric_long', 'bic_b_1', 'tric_long_1', 'tric_long_2', 'tric_long_3',
                            'tric_long_4']:
            nethumeralvector = muscleattachments[0] - muscleattachments[1]
            nethumeralvector = nethumeralvector / np.linalg.norm(nethumeralvector)

        # If no attachments on humerus, the net humeral vector is [0, 0, 0]
        else:
            nethumeralvector = np.zeros(3)

        # Each muscle's XYZ vector contributions to the GH JRF (expressed as external)
        gh_jrf_muscle[:, muscle] = nethumeralvector

    return gh_jrf_muscle, allmusclenames


def createmotfile_forces(forces, position, torques, genericmotfile, outputfilename):

    """
    Create an OpenSim compatible .mot file based on force, torque, and position inputs.

    Arguments:
        forces: Forces in Newtons (x, y, z).
        position: Positions in meters (x, y, z).
        torques: Torques in Newton-meters (x, y, z).
        genericmotfile: a template mot file.
        outputfilename: output file name.
    """

    # Using generic file as a template
    nheaders = 7
    rows = readmotfile(genericmotfile, nheaders)

    # Define external loads vector
    externalloads_vector = np.concatenate((forces, position, torques))

    # Adding time column
    time = 0

    # Save .mot file
    with open(outputfilename, 'w') as file:
        writer = csv.writer(file, delimiter='\t', lineterminator='\n')

        for row, content in enumerate(rows):
            if row == 0:
                content[0] = outputfilename
                pass

            if row > nheaders-1:
                content = np.insert(externalloads_vector, 0, time)
                time += 0.01

            writer.writerow(content)


def createmot_staticpose(kinematicsdata, genericmotfile, outputfilename):

    """
    Create an OpenSim compatible .mot file based on static pose (10 frames).

    Arguments:
        kinematicsdata: static pose.
        genericmotfile: a template mot file.
        outputfilename: output file name.
    """

    # Using generic file as a template
    nheaders = 11
    rows = readmotfile(genericmotfile, 11)

    # Adding time column
    time = 0

    # Save .mot file
    with open(outputfilename, 'w') as file:
        writer = csv.writer(file, delimiter='\t', lineterminator='\n')

        for row, content in enumerate(rows):
            if row == 0:
                content[0] = outputfilename
                pass

            if row > nheaders-1:
                content = np.insert(kinematicsdata, 0, time)
                time += 0.01

            writer.writerow(content)


def osim_muscleinfo(model, state):

    """
    Use OpenSim's API to output muscle information (maximum muscle isometric forces and moments).

    """

    nMuscles = model.getForceSet().getSize()
    nDofs = model.getCoordinateSet().getSize()
    r = np.zeros((nDofs, nMuscles))
    Fiso = np.zeros((nMuscles, nMuscles))

    # Realize the state explicitly before querying muscle quantities. equilibrateMuscles
    # expects a state realized to the Velocity stage, and moment arms require Position.
    model.realizePosition(state)
    model.realizeVelocity(state)

    for muscle in range(nMuscles):
        activemuscle = model.getMuscles().get(muscle)
        activemuscle.setActivation(state, 1)
        model.equilibrateMuscles(state)
        Fiso[muscle, muscle] = activemuscle.getTendonForce(state)
        for dof in range(nDofs):
            r[dof, muscle] = activemuscle.computeMomentArm(state, model.getCoordinateSet().get(dof))
        activemuscle.setActivation(state, 0)
        model.equilibrateMuscles(state)
    maxmoment = r @ Fiso

    return Fiso, maxmoment


def osim_run_id(model, id_setupfile, posefile, loadsfile, ti, tf, resultsdir):

    """
    Run OpenSim's Inverse Dynamics tool.

    """

    # Setup ID tool
    idTool = osim.InverseDynamicsTool(id_setupfile)
    idTool.setModel(model)
    idTool.setModelFileName(model.getName())
    idTool.setName(os.path.basename(posefile))
    idTool.setCoordinatesFileName(posefile)
    idTool.setStartTime(ti)
    idTool.setEndTime(tf)
    idTool.set_results_directory(resultsdir)
    idTool.setOutputGenForceFileName(os.path.basename(posefile)[:-8] + 'id.sto')
    idTool.setExternalLoadsFileName(resultsdir + '/' + loadsfile)

    # Save the settings in a setup file
    idsetupfile_task = os.path.basename(posefile)[:-9] + '-idsetup' + '.xml'
    idTool.printToXML(idsetupfile_task)
    shutil.move(os.path.abspath(idsetupfile_task), resultsdir + '/' + idsetupfile_task)

    # Run ID
    idTool.run()


def verify_tool_model(tool, setupfile, modelfile):

    """
    Check that an OpenSim AnalyzeTool has loaded the model we intend it to use.

    AnalyzeTool loads a model at construction from the <model_file> entry of its setup
    file, and setModelFilename() only changes the recorded name - it does not reload.
    A setup file naming a different model therefore causes the analysis to run on that
    model with no warning, which previously left Static Optimization and Joint Reaction
    Analysis running on different models. Relative paths in a setup file are resolved by
    OpenSim against the directory containing that setup file, so they are resolved the
    same way here.
    """

    declared = tool.getModelFilename()
    if not os.path.isabs(declared):
        declared = os.path.join(os.path.dirname(os.path.abspath(setupfile)), declared)
    if os.path.realpath(declared) != os.path.realpath(modelfile):
        raise RuntimeError(
            'Setup file {} names model\n  {}\nbut the analysis must run on\n  {}\n'
            'Fix <model_file> in the setup file; setModelFilename() will not reload it.'
            .format(setupfile, os.path.realpath(declared), os.path.realpath(modelfile)))


def osim_run_jra(modelfile, so_setupfile, jra_setupfile, posefile, loadsfile, ti, tf, resultsdir, forcesetfile):

    """
    Perform JRA in OpenSim without any muscles to get the JRFs due to external load + gravity alone
    This requires first running SO without any muscles, then JRA (see checks folder along with following forum post)
    https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=91&t=14425&p=0&start=0&view=&sid=c13ae20f253b57c095afe12b4b98f0be
    
    """

    # Load model, initialize, and identify relevant model parameters
    model_nomuscle = osim.Model(modelfile)
    state = model_nomuscle.initSystem()

    # Setup Analyze Tool for SO
    soTool = osim.AnalyzeTool(so_setupfile)
    verify_tool_model(soTool, so_setupfile, modelfile)
    forceset = osim.ArrayStr() # Deal with whitespace pathways
    forceset.append(os.path.relpath(forcesetfile, resultsdir))
    soTool.setForceSetFiles(forceset)  # reserve actuators
    soTool.setModelFilename(modelfile)
    soTool.setName(os.path.basename(posefile[:-9]))
    soTool.setCoordinatesFileName(posefile)
    soTool.setInitialTime(ti)
    soTool.setFinalTime(tf)
    soTool.setResultsDir(resultsdir)
    soTool.setExternalLoadsFileName(resultsdir + '/' + loadsfile)
    soTool.getAnalysisSet().get('StaticOptimization').setStartTime(ti)
    soTool.getAnalysisSet().get('StaticOptimization').setEndTime(tf)

    # Save the settings in a setup file
    so_setupfiletask = os.path.basename(posefile)[:-9] + '-sosetup' + '.xml'
    soTool.printToXML(so_setupfiletask)
    shutil.move(os.path.abspath(so_setupfiletask), resultsdir + '/' + so_setupfiletask)

    # Run SO
    soTool.run()

    # Setup Analyze Tool for JRA
    soforcefile = posefile[:-9] + '_StaticOptimization_force.sto'
    jraTool = osim.AnalyzeTool(jra_setupfile)
    verify_tool_model(jraTool, jra_setupfile, modelfile)
    jraTool.setForceSetFiles(forceset)  # reserve actuators
    jraTool.setModelFilename(modelfile)
    jraTool.setName(os.path.basename(posefile[:-9]))
    jraTool.setCoordinatesFileName(posefile)
    jraTool.setInitialTime(ti)
    jraTool.setFinalTime(tf)
    jraTool.setResultsDir(resultsdir)
    jraTool.setExternalLoadsFileName(resultsdir + '/' + loadsfile)
    jraTool.getAnalysisSet().get('JointReaction').setStartTime(ti)
    jraTool.getAnalysisSet().get('JointReaction').setEndTime(tf)
    jraTool_downcast = osim.JointReaction.safeDownCast(jraTool.getAnalysisSet().get('JointReaction'))
    jraTool_downcast.setForcesFileName(soforcefile)

    # Save the settings in a setup file
    jra_setupfiletask = os.path.basename(posefile)[:-9] + '-jrasetup' + '.xml'
    jraTool.printToXML(jra_setupfiletask)
    shutil.move(os.path.abspath(jra_setupfiletask), resultsdir + '/' + jra_setupfiletask)

    # Run JRA
    jra_setupfilepath = resultsdir + '/' + jra_setupfiletask
    jraTool = osim.AnalyzeTool(jra_setupfilepath)
    verify_tool_model(jraTool, jra_setupfilepath, modelfile)
    jraTool.run()


def readmotfile(motfile, skiprows, dataarray=None):

    """
    Read .mot data.

    Arguments:
        motfile: pathway and name of the file.
        skiprows: number of rows in headers to skip
        @:param datarray: default to return as list, otherwise return as numpy array.

    Returns:
        Data in numpy array format or in list format.

    """

    with open(motfile, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        rows = []
        for row in reader:
            rows.append(row)

    data = pd.DataFrame(rows[skiprows:])
    data = data.astype(float)
    data_np = data.to_numpy()

    if dataarray is None:
        return rows
    else:
        return data_np


def scap_rotation(model, state):

    """
    Define scapular coordinate system based on model marker locations. Output transformation matrices.

    """

    # Obtain marker landmarks in given pose
    AA = model.getMarkerSet().get('AA').getLocationInGround(state)
    TS = model.getMarkerSet().get('TS').getLocationInGround(state)
    AI = model.getMarkerSet().get('AI').getLocationInGround(state)
    GH = model.getMarkerSet().get('GH').getLocationInGround(state)
    Glenoid = model.getMarkerSet().get('Glenoid').getLocationInGround(state)
    AA = np.array([AA[0], AA[1], AA[2]])
    TS = np.array([TS[0], TS[1], TS[2]])
    AI = np.array([AI[0], AI[1], AI[2]])
    GH = np.array([GH[0], GH[1], GH[2]])
    Glenoid = np.array([Glenoid[0], Glenoid[1], Glenoid[2]])
    origin = Glenoid

    # Define scapular coordinate system
    # X+: Lateral;  Y+: Superior;  Z+: Backwards
    xcs = (Glenoid - TS) / np.linalg.norm(Glenoid - TS)
    temp = (AI - TS) / np.linalg.norm(AI - TS)
    zcs = np.cross(temp, xcs) / np.linalg.norm(np.cross(temp, xcs))
    ycs = np.cross(zcs, xcs)
    R = np.column_stack((xcs, ycs, zcs))
    Rinv = np.linalg.inv(R)

    return R, Rinv, origin


def visualizeoutputs(moments, dofnames, maxforces, musclenames, maxmoments, gh_jrfs, savepath):

    """
    Visualize all data: inverse dynamics moments, maximum muscle forces, maximum muscle moments, GH joint reaction forces.

    """

    # Plot ID moments
    plt.barh(range(len(moments[6:])), moments[6:], tick_label=dofnames[6:])
    plt.axvline(x=0, color='black', linestyle='--')
    plt.xlabel('Moment (Nm)')
    plt.ylabel('DOF')
    plt.tight_layout()
    plt.savefig(savepath + '-loads.png', format='png')
    plt.close()

    # Plot max muscle forces
    plt.barh(range(len(maxforces)), maxforces.diagonal(), tick_label=musclenames)
    plt.xlabel('Max Isometric Force (N)')
    plt.ylabel('Muscle')
    plt.yticks(fontsize=6)
    plt.tight_layout()
    plt.savefig(savepath + '-Fiso.png', format='png')
    plt.close()

    # Plot max muscle moments
    for i in range(6, 17):
        plt.barh(range(len(maxmoments[i, :])), maxmoments[i, :], tick_label=musclenames)
        plt.axvline(x=0, color='black', linestyle='--')
        plt.xlabel('Max Moment (Nm)')
        plt.ylabel('Muscle')
        plt.yticks(fontsize=6)
        plt.title(dofnames[i])
        plt.tight_layout()
        plt.savefig(savepath + '-maxmoments-' + dofnames[i] + '.png', format='png')
        plt.close()

    # Plot stability ratios
    stabilityratio = ['Sy', 'Sz']
    for i in range(2):
        stability = gh_jrfs[i+1, :] / abs(gh_jrfs[0, :])
        stability = np.nan_to_num(stability, nan=0)
        plt.barh(range(len(stability)), stability, tick_label=musclenames)
        plt.axvline(x=0, color='black', linestyle='--')
        plt.xlabel('Stability Ratio')
        plt.ylabel('Muscle')
        plt.yticks(fontsize=6)
        plt.title(stabilityratio[i])
        plt.tight_layout()
        plt.savefig(savepath + '-stability-' + stabilityratio[i] + '.png', format='png')
        plt.close()