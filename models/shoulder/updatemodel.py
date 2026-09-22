import numpy as np
import opensim as osim


def updatetothelen(modelfilepath):
    """
    Converts .osim model with Schutte 1993 muscle model to Thelen 2003 muscle model.

    :param modelfilepath: .osim file path
    """

    # Define the word to find and replace
    word_to_find = 'Schutte1993Muscle_Deprecated'
    replacement_word = 'Thelen2003Muscle'
    words_to_check = ['Parameter used in time constant of ramping up', 'activation1', 'activation2']

    # New filename
    modelfilepath_new = modelfilepath[:-5] + '-Thelen' + modelfilepath[-5:]

    # Open the file in read mode
    with open(modelfilepath, 'r') as file:
        content = file.read()

    # Replace all instances of the word
    modified_content = content.replace(word_to_find, replacement_word)

    # Open the file again in write mode to save the modified content
    with open(modelfilepath_new, 'w') as file:
        file.write(modified_content)

    # Read the content of the file
    with open(modelfilepath_new, 'r') as file:
        lines = file.readlines()

    # Filter out lines that contain any of the words in words to check
    filtered_lines = [line for line in lines if not any(word in line for word in words_to_check)]

    # Open the output file in write mode to save the filtered lines
    with open(modelfilepath_new, 'w') as output_file:
        output_file.writelines(filtered_lines)


def simplifydsem(modelfilepath):
    """
    Converts DSEM to simplified model version with maximum 3 elements per muscle.
    Strength of muscles removed is added to muscles kept.

    :param modelfilepath: .osim file path
    """

    # Load and initialize model
    model = osim.Model(modelfilepath)
    nMuscles = model.getMuscles().getSize()

    # List of muscles to keep
    keepmuscles = ['conoid_lig', 'trap_scap_3', 'trap_scap_7', 'trap_scap10', 'trap_clav_1', 'lev_scap_1', 'pect_min_2',
                   'rhomboid_2', 'rhomboid_4', 'serr_ant_2', 'serr_ant_6', 'serr_ant10', 'delt_scap_2', 'delt_scap10',
                   'delt_clav_2', 'coracobr_2', 'infra_2', 'infra_4', 'ter_min_2', 'ter_maj_3', 'supra_1', 'supra_3',
                   'subscap_2', 'subscap_5', 'subscap_9', 'bic_l', 'bic_b_1', 'tric_long_2', 'lat_dorsi_2',
                   'lat_dorsi_4', 'lat_dorsi_6', 'pect_maj_t_3', 'pect_maj_t_5', 'pect_maj_c_2', 'tric_med_4',
                   'brachialis_4', 'brachiorad_2', 'pron_teres_1', 'pron_teres_2', 'supinator_3', 'pron_quad_2',
                   'tric_lat_4', 'anconeus_2']
    revisedmusclenames = ['conoid_lig', 'trap_scap_S', 'trap_scap_M', 'trap_scap_I', 'trap_clav', 'lev_scap',
                          'pect_min', 'rhomboid_S', 'rhomboid_I', 'serr_ant_I', 'serr_ant_M', 'serr_ant_S',
                          'delt_scap_P', 'delt_scap_M', 'delt_clav_A', 'coracobr', 'infra_I', 'infra_S', 'ter_min',
                          'ter_maj', 'supra_P', 'supra_A', 'subscap_S', 'subscap_M', 'subscap_I', 'bic_l', 'bic_b',
                          'tric_long', 'lat_dorsi_S', 'lat_dorsi_M', 'lat_dorsi_I', 'pect_maj_t_I', 'pect_maj_t_M',
                          'pect_maj_c_S', 'tric_med', 'brachialis', 'brachiorad', 'pron_teres_H', 'pron_teres_U',
                          'supinator', 'pron_quad', 'tric_lat', 'anconeus']

    # List of muscles at the end of each muscle group (required for summing up forces)
    edgemuscles = ['conoid_lig', 'trap_scap_5', 'trap_scap_8', 'trap_scap11', 'trap_clav_2', 'lev_scap_2', 'pect_min_4',
                   'rhomboid_3', 'rhomboid_5', 'serr_ant_4', 'serr_ant_8', 'serr_ant12', 'delt_scap_3', 'delt_scap11',
                   'delt_clav_4', 'coracobr_3', 'infra_3', 'infra_6', 'ter_min_3', 'ter_maj_4', 'supra_2', 'supra_4',
                   'subscap_3', 'subscap_6', 'subscap11', 'bic_l', 'bic_b_1', 'tric_long_4', 'lat_dorsi_2',
                   'lat_dorsi_4', 'lat_dorsi_6', 'pect_maj_t_4', 'pect_maj_t_6', 'pect_maj_c_2', 'tric_med_5',
                   'brachialis_7', 'brachiorad_3', 'pron_teres_1', 'pron_teres_2', 'supinator_5', 'pron_quad_3',
                   'tric_lat_5', 'anconeus_5']

    # List of all muscle names and forces
    allmuscles = []
    allF0 = []
    edgemuscles_index = []
    for muscle in range(nMuscles):
        # Get muscle information
        activemuscle = model.getMuscles().get(muscle)
        activemuscle_name = activemuscle.getName()
        activemuscle_F0 = activemuscle.getMaxIsometricForce()

        # Append list
        allmuscles.append(activemuscle_name)
        allF0.append(activemuscle_F0)

        if activemuscle_name in edgemuscles:
            edgemuscles_index.append(allmuscles.index(activemuscle_name))

    # Updating force of muscles to be kept
    counter = 0
    for muscle in allmuscles:
        activemuscle = model.getMuscles().get(muscle)

        if muscle in keepmuscles:
            if counter == 0:
                # First group runs from the start of the list to its edge element.
                # (The old slice was empty and only happened to give the right answer
                # because the first kept element, conoid_lig, has F0 = 0.)
                updforce = np.sum(allF0[0:edgemuscles_index[counter] + 1])
            else:
                updforce = np.sum(allF0[edgemuscles_index[counter - 1] + 1:edgemuscles_index[counter] + 1])

            activemuscle.setMaxIsometricForce(updforce)
            counter += 1

        else:
            pass

    # Removing muscles
    for muscle in allmuscles:
        activemuscle = model.getMuscles().get(muscle)

        if muscle in keepmuscles:
            keepmuscles_index = keepmuscles.index(muscle)
            activemuscle.setName(revisedmusclenames[keepmuscles_index])

        else:
            removemuscle_index = model.getMuscles().getIndex(muscle)
            model.updForceSet().remove(removemuscle_index)

    # Saving model
    model.finalizeConnections()
    modelfile_new = modelfilepath[:-5] + '-simplified' + modelfilepath[-5:]
    model.printToXML(modelfile_new)


if __name__ == "__main__":

    modelfile = 'das3_clav_scap_orig_PhD.osim'
    updatetothelen(modelfile)
    simplifydsem(modelfile[:-5] + '-Thelen' + modelfile[-5:])
