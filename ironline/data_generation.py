import os
import time
import subprocess
import numpy as np
# from typing import List, Dict


def parameter_tensor(spin_list: list, defpar_list: list, inc_list: list, ode_list: list) -> np.ndarray:
    """Generate a tensor of parameters to be fed into the main.cpp file

    Args:
        spin_list (list): List of spins to be used
        defpar_list (list): List of defpar values to be used
        inc_list (list): List of inclination angles to be used

    Returns:
        np.ndarray: Tensor of parameters to be fed into the main.cpp file
    """

    # tensor = np.zeros((len(spin_list), len(defpar_list), len(inc_list), 3))
    # for i, spin in enumerate(spin_list):
    #     for j, defpar in enumerate(defpar_list):
    #         for k, angle in enumerate(inc_list):
    #             tensor[i, j, k] = [spin, defpar, angle]
    
    # Generate 4D tensor of parameters to be fed into the main.cpp file
    tensor = np.zeros((len(spin_list), len(defpar_list), len(inc_list), len(ode_list), 4))
    for i, spin in enumerate(spin_list):
        for j, defpar in enumerate(defpar_list):
            for k, angle in enumerate(inc_list):
                for l, ode in enumerate(ode_list):
                    tensor[i, j, k, l] = [spin, defpar, angle, ode]
    return tensor


def command_generation(params: dict, executable_path: str,
                       errtol: float, rstep: float, pstep: float,
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
    defpar = params['defpar']
    angle = params['inc']
    ode_solver = params['ode_solver']

    # Create c++ command to be executed
    command = [executable_path, str(spin), str(defpar), str(angle), str(
        errtol), str(rstep), str(pstep), str(progress_check), str(ode_solver)]
    command = ' '.join(command)

    return command


def data_generation(spin_list: list, defpar_list: list, inc_list: list,
                    executable_path: str,
                    errtol: float = 1.0e-8, rstep: float = 1.01, pstep: float = 0.007853981634, progress_check: int = 0, ode_solver: int = 0) -> None:
    """Generate data using the main.cpp file

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
        ode_solver (int, optional): ODE solver to be used. 0 for RK45, 1 for RKN. Defaults to 0 (RK45).
    """

    # Generate tensor of parameters to be fed into the main.cpp file
    params = parameter_tensor(spin_list, defpar_list, inc_list, [ode_solver])

    # Start timer
    start = time.time()
    mid = start

    # Loop through the tensor of parameters and generate data
    for i in range(params.shape[0]):
        for j in range(params.shape[1]):
            for k in range(params.shape[2]):
                param_dict = {
                    'spin': params[i, j, k, 0], 'defpar': params[i, j, k, 1], 'inc': params[i, j, k, 2]}
                command = command_generation(
                    param_dict, executable_path, errtol, rstep, pstep, progress_check, ode_solver)
                os.system(command)

                # Estimate remaining time
                time_taken = time.time() - mid
                mid = time.time()
                total_iter = params.shape[0]*params.shape[1]*params.shape[2]
                current_iter = i * \
                    params.shape[1]*params.shape[2] + j*params.shape[2] + k + 1
                remaining = (total_iter - current_iter)*time_taken

                print(
                    f'Current iteration: {current_iter}/{total_iter}, {current_iter/total_iter*100} percent finished')
                print(
                    f'Estimated time remaining: {remaining/60} minutes or {remaining/3600} hours')

                # if output_path is None:
                #     output_path = os.getcwd()
                # # Move data to output directory
                # filename = 'iron_a{:.3f}.def{:.2f}.i{:.2f}.dat'.format(
                #     param_dict['spin'], param_dict['defpar'], param_dict['inc'])
                # os.system('mv {} {}'.format(filename, output_path/filename))


def data_generation_para(spin_list: list, defpar_list: list, inc_list: list, ode_solver: list,
                         executable_path: str, errtol: float = 1.0e-8, rstep: float = 1.01, pstep: float = 0.007853981634, progress_check: int = 0) -> None:
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
    params = parameter_tensor(spin_list, defpar_list, inc_list, ode_solver)

    # Spawn processes each running a different set of parameters
    for i in range(params.shape[0]):
        for j in range(params.shape[1]):
            for k in range(params.shape[2]):
                for l in range(params.shape[3]):
                    param_dict = {
                        'spin': params[i, j, k, l, 0], 'defpar': params[i, j, k, l, 1], 'inc': params[i, j, k, l, 2], 'ode_solver': params[i, j, k, l, 3]}

                    command = command_generation(
                        param_dict, executable_path, errtol, rstep, pstep, progress_check)
                    print(command)
                    process_list.append(subprocess.Popen(command, shell=True))

    # Wait for processes to finish
    exit_codes = [p.wait() for p in process_list]
    print(exit_codes)


if __name__ == '__main__':
    spins = [0.10, 0.50, 0.998]
    defpars = [0.0, 5.00, 10.00]
    incs = [45.0]

    # spins = [0.10, 0.50, 0.998]
    # defpars = [5.00, 10.00]
    # incs = [20.0, 45.0, 70.0]

    odes = [1]

    PATH_TO_EXECUTABLE = r'/workspaces/raytransfer'
    PATH_TO_OUTPUT = r'C:\Users\WalkerXin\Documents\Scripts\raytransfer\ironline\data'
    os.chdir(PATH_TO_EXECUTABLE)
    print(os.getcwd())

    # Estimate time taken, with each run taking ~ 8 minutes
    size = len(spins)*len(defpars)*len(incs)*len(odes)
    estimate = 15
    print(f'Estimated time taken: {size*estimate} minutes or {size*estimate/60} hours')

    # Prompt user to proceed, default: y
    check = input('Proceed? (y/n): ') or 'y'

    if check == 'Y' or 'y':
        pass
    else:
        exit()

    # Generate data
    for defpar in defpars:
            data_generation_para(
                spins, [defpar], incs, [0], r'ironline/main',
                errtol=1.0e-10, rstep=1.008, pstep=2*np.pi/800, progress_check=2)

    # for ode in odes:
    #         if ode == 0:
    #             for defpar in defpars:
    #                 data_generation_para(
    #                     spins, [defpar], incs, [ode], r'/workspaces/raytransfer/ironline/main',
    #                     errtol=1.0e-10, rstep=1.008, pstep=2*np.pi/800, progress_check=0)
    #         else:
    #             for defpar in defpars:
    #                 data_generation_para(
    #                     spins, [defpar], incs, [ode], r'C:\Users\WalkerXin\Documents\Scripts\raytransfer\ironline\main.exe',
    #                     errtol=1.0e-10, rstep=1.008, pstep=2*np.pi/800, progress_check=2)

    

    # Turn off computer 60 seconds after completion
    # os.system('shutdown /s /t 60')

    # # Ask for user input to cancel shutdown
    # check = input('Shutdown in 60 seconds. Cancel? (y/n): ') or 'n'
    # if check == 'Y' or 'y':
    #     os.system('shutdown /a')
    # else:
    #     pass
