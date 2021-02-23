# This file contains functions that facilitate parallel programming

import math
from solver_functions import *

####### Defining MPI functions ############################################################################
###################################################################################################
###################################################################################################


def divide_array(data, p_n, p_i):
    ### gets and array and calculates share each processor
    data_n = len(data)
    share_n = math.ceil(data_n / p_n);
    data_p_s = p_i*share_n
    data_p_e = min((p_i+1)*share_n, data_n)

    return data_p_s, data_p_e



def collect_array(comm, data_2, p_n, p_i):
    # gets diffrents parts of an array from all processors and appends and returns the total array
    if p_i != 0:
        comm.send(data_2, dest=0)
    else:
        for i in range(1,p_n):
            data_2 = data_2 + comm.recv(source=i)

    if p_n > 1:
        data_2 = comm.bcast(data_2, root=0)

    return data_2






def collect_maximums(comm, p_n, p_i, max_f, max_solution, max_combination, max_feasible_interval, max_interval):
    # this function gets array of intervals with their maximum values from different processors and
    # finds one array that contains intervals with the overall maximum values

    if p_i != 0: # except for process zero send arrays to process zero
        comm.send(max_f, dest=0)
        comm.send(max_solution, dest=0)
        comm.send(max_combination, dest=0)
        comm.send(max_feasible_interval, dest=0)
        comm.send(max_interval, dest=0)

    else:
        for i in range(1,p_n): # process 0 gets all the maximum intervals and updates its maximum interval
            max3_f = comm.recv(source=i)
            max3_solution = comm.recv(source=i)
            max3_combination = comm.recv(source=i)
            max3_feasible_interval = comm.recv(source=i)
            max3_interval = comm.recv(source=i)
            # updating intervals using intervals from other processors
            max_f, max_solution, max_combination, max_feasible_interval, max_interval = maximum_2(max_f, max_solution, max_combination, max_feasible_interval, max_interval, max3_f, max3_solution, max3_combination, max3_feasible_interval, max3_interval)

    return max_f, max_solution, max_combination, max_feasible_interval, max_interval
