import os
import subprocess
import numpy as np
# from typing import List, Dict


def parameter_tensor(spin_list: list, inc_list: list, defpar_list: list) -> np.ndarray:
    """Generate a tensor of parameters to be fed into the main.cpp file

    Args:
        spin_list (list): List of spins to be used
        inc_list (list): List of inc values to be used
        defpar_list (list): List of defpar values to be used

    Returns:
        np.ndarray: Tensor of parameters to be fed into the main.cpp file
    """

    tensor = np.zeros((len(spin_list), len(inc_list), len(defpar_list), 3))
    for i, spin in enumerate(spin_list):
        for j, inc in enumerate(inc_list):
            for k, angle in enumerate(defpar_list):
                tensor[i, j, k] = [spin, inc, angle]

    return tensor


def command_generation(params: dict, executable_path: str,
                       gerrtol: float, rerrtol: float,
                       progress_check: int) -> str:
    """Generate data using the main.cpp file

    Args:
        params (dict): Dictionary of parameters to be fed into the main.cpp file
        executable_path (str): Path to the main.cpp executable
        errtol (float): Error tolerance
        rstep (float): Step size for robs
        pstep (float): Step size for pobs
        progress_check (int): Number of iterations to check progress, zero means no progress check

    Returns:
        str: c++ command to be executed
    """

    # Retrieve parameters to be passed as argv to the main.cpp file
    spin = params['spin']
    inc = params['inc']
    defpar = params['defpar']

    # Create c++ command to be executed
    command = [executable_path, str(spin), str(inc), str(
        defpar), str(gerrtol), str(rerrtol), str(progress_check)]
    command = ' '.join(command)

    return command


def data_generation(spin_list: list, inc_list: list, defpar_list: list,
                         executable_path: str, gerrtol: float = 1.0e-8, rerrtol: float = 1.0e-8, progress_check: int = 0) -> None:
    """Generate data using the main.cpp file with parallel processing

    Args:
        spin_list (list): List of spins to be used
        defpar_list (list): List of defpar values to be used
        inc_list (list): List of inclination angles to be used
        executable_path (str): Path to the main.cpp executable
        output_path (str): Path to the output directory. Defaults to None.
        errtol (float, optional): Error tolerance. Defaults to 1.0e-8.
        rstep (float, optional): Step size for robs. Defaults to 1.01.
        pstep (float, optional): Step size for pobs. Defaults to 2*Pi/800.
        progress_check (int, optional): Number of iterations to check progress. Defaults to 0 (no progress check).
    """
    process_list = []

    # Generate tensor of parameters to be fed into the main.cpp file
    params = parameter_tensor(spin_list, inc_list, defpar_list)

    # Spawn processes each running a different set of parameters
    for i in range(params.shape[0]):
        for j in range(params.shape[1]):
            for k in range(params.shape[2]):
                param_dict = {
                    'spin': params[i, j, k, 0], 'inc': params[i, j, k, 1], 'defpar': params[i, j, k, 2]}

                command = command_generation(
                    param_dict, executable_path, gerrtol, rerrtol, progress_check)
                # print(command)
                process_list.append(subprocess.Popen(command))

    # Wait for processes to finish
    exit_codes = [p.wait() for p in process_list]
    print(exit_codes)
    
    
def data_generation_tensor(input_tensor, executable_path: str, gerrtol: float = 1.0e-8, rerrtol: float = 1.0e-8, progress_check: int = 0, size_check = None):
    if size_check is None:
        size_check = input_tensor.shape[0]
    
    # Check size of input tensor
    if input_tensor.shape[0] != size_check:
        raise ValueError('Input tensor does not match size check')
    
    process_list = []
    
    # Spawn processes each running a different set of parameters
    for element in input_tensor:
                param_dict = {
                    'spin': element[0], 'inc': element[1], 'defpar': element[2]}

                command = command_generation(
                    param_dict, executable_path, gerrtol, rerrtol, progress_check)
                # print(command)
                process_list.append(subprocess.Popen(command))
                
    # Wait for processes to finish
    exit_codes = [p.wait() for p in process_list]
    print(exit_codes)
        


if __name__ == '__main__':
    spins = [0.5, 0.9, 0.998]
    incs = [45, 70]
    defpars = [0.0, 1.0, 5.0, 10.0]
    dissect = 8
    # defpars = [0.0]
    
    PATH_TO_EXECUTABLE = r'C:\Users\WalkerXin\Documents\Scripts\raytransfer'
    PATH_TO_OUTPUT = r'C:\Users\WalkerXin\Documents\Scripts\raytransfer\photons'
    os.chdir(PATH_TO_EXECUTABLE)
    
    # Make sure that photons directory exists
    if not os.path.exists(PATH_TO_OUTPUT):
        os.makedirs(PATH_TO_OUTPUT)

    # Estimate time taken, with each run taking ~ 8 minutes
    size = len(spins)*len(incs)*len(defpars)
    estimate = 120
    print(
        f'Estimated time taken: {size*estimate} minutes or {size*estimate/60} hours')

    # Prompt user to proceed, default: y
    check = input('Proceed? (y/n): ') or 'y'

    if check == 'Y' or 'y':
        pass
    else:
        exit()
    
    # Empty array that stores 3-tuple
    space = np.array(np.zeros((len(spins), len(incs), len(defpars))), dtype=object)
    for i, spin in enumerate(spins):
        for j, inc in enumerate(incs):
            for k, angle in enumerate(defpars):
                space[i, j, k] = (spin, inc, angle)
    # Flatten array
    space = space.flatten()

    # Split array into chunks of size dissect
    chunks = []
    for i in range(0, len(space), dissect):
        chunk = space[i:i+dissect]
        chunks.append(chunk)
    
    for chunk in chunks:
        try:
            data_generation_tensor(chunk, r'C:\Users\WalkerXin\Documents\Scripts\raytransfer\main.exe',
                               gerrtol=1.0e-6, rerrtol=1.0e-7, progress_check=1, size_check=dissect)
        except ValueError:
            print('Size check failed')
            continue
    
    # data_generation_para(spins, incs, defpars, r'C:\Users\WalkerXin\Documents\Scripts\raytransfer\main.exe',
    #                           gerrtol=1.0e-6, rerrtol=1.0e-7, progress_check=1)

    # Turn off computer 60 seconds after completion
    # os.system('shutdown /s /t 60')

    # # Ask for user input to cancel shutdown
    # check = input('Shutdown in 60 seconds. Cancel? (y/n): ') or 'n'
    # if check == 'Y' or 'y':
    #     os.system('shutdown /a')
    # else:
    #     pass
