import os
def evaluate_river_networks(sFilename_conceptual_river_network_in, sFilnanme_actual_river_network_in, sFilename_output_in=None):
    """
    Evaluate the river network by comparing the conceptual river network with the actual river network.
    Args:
        sFilename_conceptual_river_network_in (str): The conceptual river network file.
        sFilnanme_actual_river_network_in (str): The actual river network file.
    """

    #check if the files exist

    if not os.path.exists(sFilename_conceptual_river_network_in):
        print(f"The file {sFilename_conceptual_river_network_in} does not exist.")
    if not os.path.exists(sFilnanme_actual_river_network_in):
        print(f"The file {sFilnanme_actual_river_network_in} does not exist.")
        return


    #because the actual river network is main channel, we should only use the main channel in the conceptual river network
    
    return