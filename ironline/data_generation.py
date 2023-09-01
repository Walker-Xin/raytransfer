# Helper function to run main.cpp across a matrix of parameters

import os
import numpy as np
import typing
import time


def parameter_tensor(spin_list: list, defpar_list: list, inc_list: list) -> np.ndarray:
    """Generate a tensor of parameters to be fed into the main.cpp file

    Args:
        spin_list (list): List of spins to be used
        defpar_list (list): List of defpar values to be used
        inc_list (list): List of inclination angles to be used

    Returns:
        np.ndarray: Tensor of parameters to be fed into the main.cpp file
    """
    tensor = np.zeros((len(spin_list), len(defpar_list), len(inc_list), 3))
    for i, spin in enumerate(spin_list):
        for j, defpar in enumerate(defpar_list):
            for k, inc in enumerate(inc_list):
                tensor[i, j, k] = [spin, defpar, inc]
    return tensor


def command_generation(params: dict, executable_path: str, errtol: float, rstep: float, pstep: float) -> str:
    """Generate data using the main.cpp file

    Args:
        params (dict): Dictionary of parameters to be fed into the main.cpp file
        output_dir (str): Directory to store the output files
        executable_path (str): Path to the main.cpp executable
    """

    # Retrieve parameters to be passed as argv to the main.cpp file
    spin = params['spin']
    defpar = params['defpar']
    inc = params['inc']

    # Create c++ command to be executed
    command = [executable_path, str(spin), str(defpar), str(inc), str(errtol), str(rstep), str(pstep)]
    command = ' '.join(command)

    return command


def data_generation(spin_list: list, defpar_list: list, inc_list: list, executable_path, errtol: float = 1.0e-8, rstep: float = 1.01, pstep: float = 0.007853981634) -> None:
    """Generate data using the main.cpp file

    Args:
        spin_list (list): List of spins to be used
        defpar_list (list): List of defpar values to be used
        inc_list (list): List of inclination angles to be used
        executable_path (str): Path to the main.cpp executable
        errtol (float, optional): Error tolerance. Defaults to 1.0e-8.
        rstep (float, optional): Step size for robs. Defaults to 1.01.
        pstep (float, optional): Step size for pobs. Defaults to 0.007853981634.
    """

    # Generate tensor of parameters to be fed into the main.cpp file
    params = parameter_tensor(spin_list, defpar_list, inc_list)

    # Start timer
    start = time.time()

    # Loop through the tensor of parameters and generate data
    for i in range(params.shape[0]):
        for j in range(params.shape[1]):
            for k in range(params.shape[2]):
                param_dict = {
                    'spin': params[i, j, k, 0], 'defpar': params[i, j, k, 1], 'inc': params[i, j, k, 2]}
                command = command_generation(
                    param_dict, executable_path, errtol, rstep, pstep)
                os.system(command)

    # End timer and print time taken
    end = time.time()
    print('Time taken: {} minutes'.format((end-start)/60))


if __name__ == '__main__':
    spin_list = [-0.50, 0.50, 0.98]
    defpar_list = [0.00, 0.50, 1.00, 2.00]
    inc_list = [30.0, 45.0, 65.0]

    # spin_list = [0.9]
    # defpar_list = [1.00]
    # inc_list = [45.0]

    path_to_executable = r'C:\Users\WalkerXin\Documents\Scripts\raytransfer\ironline'
    os.chdir(path_to_executable)

    # Estimate time taken, with each run taking ~ 8 minutes
    size = len(spin_list)*len(defpar_list)*len(inc_list)
    print('Estimated time taken: {} minutes or {} hours'.format(size*8, size*8/60))

    # Prompt user to proceed
    check = input('Proceed? (y/n): ')
    if check == 'Y' or 'y':
        pass
    else:
        exit()

    # Generate data
    data_generation(spin_list, defpar_list, inc_list, 'main.exe', errtol=1.0e-8, rstep=1.01, pstep=2*np.pi/800)
    pass
