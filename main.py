#### Analytical Optimizer main file
# to run this code you need these modules installed:
# Sympy
# pickle
# tabulate

# for python: python -m pip install antlr4-python3-runtime # for parse_latex() command
# or for conda: conda install -c carta antlr4-python3-runtime  # for parse_latex() command
# if you get A"NTLR runtime and generated code versions disagree: 4.8!=4.7.1" use "python -m pip install antlr4-python3-runtime==4.7.1"

# mpi4py (only if you have MPI setup and want to run in parallel)

# the program is designed to run the calculations in parallel if run by MPI but it can also be run normally.
#### To run in parallel:
# First you need to have MPI setup on your computer
# go to commandline and navigate to the code folder (no need if you are using command line in your editor)
# Run the command below:
# mpirun.openmpi -n 4 python3 main.py
# -n determines the number processors to be used. On a single computer with a cpu with 4 cores "-n 4" is the fastest.

#fgdfg

from MPI_functions import *

from sympy.parsing.latex import parse_latex
from sympy import  *
import itertools
import pickle
from tabulate import tabulate


from sympy import symbols, IndexedBase
o,s,m = symbols('o s m', real=True)   # constants
Q = IndexedBase('Q', real=True)  #main indexed symbol



###### To load MPI if it exists ###################
try:
    # try to load MPI and run in parallel
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    p_i = comm.Get_rank()
    p_n = comm.Get_size()
except:
    comm = None
    p_i = 0
    p_n = 1


####### making output folder ############################################################################
from pathlib import Path
Path("./output").mkdir(parents=True, exist_ok=True)



####### reading the function to be optimized ############################################################################
###################################################################################################
# function and constraint are defined in "input.py" and loaded here


exec(open("./run_options.py").read())

if p_i == 0:  # Only process 0 analyses the input and then broadcast variable so that every one has consistent translation between original variables abd Working variabale Q[i]
    if load_input: # load input from the output of the previous run
        with open('./output/1_input.txt', 'rb') as fp:
            f, c_exp, c_exp_eq, Constants, c2_exp, c2_exp_eq, var_list = pickle.load(fp)
    else:
        ####### reading and converting equations to sympy (indexed variables) ############################################################################
        ###################################################################################################

        # set Function, Constraints and Constants_Constraints using input.py
        exec(open("./input.py").read())

        #converting main function f to sympy format,
        #replacing varibales with indexed ones and making array of
        #variables to make iterating over variables easy

        f = parse_latex(Function) #conversion from latec to sympy expression
        var_list = list(f.free_symbols) #finding the independent variables used

        #converting constraints to sympy format and substituting variables by indexed ones
        const_n = len(Constraints)
        c_exp=[0 for i in range(const_n)] #empty list to store constraints
        for const_i in range(const_n):
            c_exp[const_i] = parse_latex(Constraints[const_i])

        c_exp_eq=[0 for i in range(const_n)] #empty function to store constraints
        for const_i in range(const_n):
            c_exp_eq[const_i]=c_exp[const_i].lhs

        #converting constraints on constants to sympy format
        const2_n = len(Constants_Constraints)
        c2_exp=[0 for i in range(const2_n)] #empty array to store constraints
        for const_i in range(const2_n):
            c2_exp[const_i] = parse_latex(Constants_Constraints[const_i])
        c2_exp_eq=[0 for i in range(const2_n)] #empty array to store constraints
        for const_i in range(const2_n):
            c2_exp_eq[const_i]=c_exp[const_i].lhs

        if p_n == 1: # saving for troubleshooting
            with open('./output/1_input.txt', 'wb') as fp:
                pickle.dump([f, c_exp, c_exp_eq, Constants, c2_exp, c2_exp_eq, var_list],fp)





# Broadcast result of reading input.py to other processors
if p_n > 1:
    if p_i != 0:
        f= 0
    f = comm.bcast(f, root=0)
    if p_i != 0:
        c_exp = 0
    c_exp = comm.bcast(c_exp, root=0)
    if p_i != 0:
        c_exp_eq = 0
    c_exp_eq = comm.bcast(c_exp_eq, root=0)
    if p_i != 0:
        Constants = 0
    Constants = comm.bcast(Constants, root=0)
    if p_i != 0:
        c2_exp = 0
    c2_exp = comm.bcast(c2_exp, root=0)
    if p_i != 0:
        c2_exp_eq = 0
    c2_exp_eq = comm.bcast(c2_exp_eq, root=0)
    if p_i != 0:
        var_list = 0
    var_list = comm.bcast(var_list, root=0)


if p_i == 0:
    print("Reading input done!")

####### Finding all possible combinations of constarints ############################################################################
###################################################################################################

constraint_combinations_all = []

if not load_constraint_combinations:
    # makes a list of all possible combinations of boundaries.
    # each combination can have 0 to var_n members.
    constraint_combinations_all = []
    for i in range(len(var_list)+1):
        constraint_combinations_all = constraint_combinations_all + list(itertools.combinations(c_exp_eq,i))

    if p_n==1 and p_i==0: # saving for troubleshooting
        with open("./output/2_constraint_combinations.txt", "wb") as fp:  # Pickling
            pickle.dump(constraint_combinations_all, fp)
else:
    with open("./output/2_constraint_combinations.txt", "rb") as fp:  # Pickling
        constraint_combinations_all = pickle.load(fp)



# for each process loading only relevant ones based on modulo
constraint_combinations = []
for i in range(len(constraint_combinations_all)):
    if i % p_n == p_i:
        constraint_combinations = constraint_combinations + [constraint_combinations_all[i]]
if p_i==0:
    print("Process", p_i, ": listing constraint combinations done!")
    print("Number of all combinations:", len(constraint_combinations_all))

# free memory
constraint_combinations_all = None


####### Finding intesect of each combination with the main function ############################################################################
###################################################################################################
###################################################################################################

if not load_intersect: # load intersections from output of the previous run

    # function intersect() from "solver_functions.py" caculates intersection and outputs lists below:
    # intersect_f: list is to save intersection of f with set of constraints
    # intersect_sol: list to save intersection of constraints with each other
    # intersect_comb_id: to keep id of combination for backward check
    intersect_f, intersect_solution, intersect_combination = intersect(f, constraint_combinations, var_list, p_i)
    # free memory
    constraint_combinations = None

    if p_n == 1: # only save in non-parallel run (saving for troubleshooting)
        with open('./output/3_intersect.txt', 'wb') as fp:
            pickle.dump([intersect_f, intersect_solution, intersect_combination], fp)
else:
    with open('./output/3_intersect.txt', 'rb') as fp:
        intersect_f, intersect_solution, intersect_combination = pickle.load(fp)



# printing all intersection results

if p_i == 0 and print_intersects:
    print("")
    print("List of intersections:")
    intersect_i = range(len(intersect_f))
    headers = ["intersect_i", "intersect_f", "intersect_sol", "intersect_comb_id"]
    table = zip(intersect_i, intersect_f, intersect_sol, intersect_comb_id)
    print(tabulate(table, headers, tablefmt="fancy_grid"))
if p_i == 0:
    print("")
    print("")




####### Finding list of optimum points of the intersects ############################################################################
###################################################################################################
###################################################################################################

# load from output of the previous run or calculate
if not load_optimum: # load optimums from output of the previous run

    # This function takes the partial derivativs of the all intersections of function f with combinations of constraints
    # with respect to all the variables usinf diff() function of Sympy then it solves the system of equations of derivatives
    # for zeros. It outputs a solution list (values of variables at optimum): "sol_list"
    # and corresponding optimum function values list: "f_opt".
    # Outputs:
    # opt_f: optimum value of f intersection with constraint combination (optimum values withing the volume of and at the boundaries)
    # opt_sol: to keep values of variables for each optimum
    # opt_comb_id: to keep combination id for backward check
    # opt_intersect_id: to keep intersection id for backward check
    # The main limitaion of the code is the ability of nonlinsolve() of Sympy
    # to sovle system of nonlinear equations. If the main function is quadratic this sytem is linear and there is certainly
    # no problem. For nonlinear functions more complicated than quadratic it may get stuk here and run forever.
    opt_f, opt_solution, opt_combination = optimum(intersect_f, intersect_solution, intersect_combination, var_list, p_n, p_i)
    #free memory
    intersect_f = None
    intersect_solution = None
    intersect_combination = None

    if p_n==1 and p_i==0:  # only save in non-parallel run (saving for troubleshooting)
        with open('./output/4_optimum.txt', 'wb') as fp:
            pickle.dump([opt_f, opt_solution, opt_combination], fp)
else:
    with open('./output/4_optimum.txt', 'rb') as fp:
        opt_f, opt_solution, opt_combination = pickle.load(fp)

# printing all optimimus
if p_i == 0 and print_optimum:
    print("")
    print("List of optimums:")
    opt_i = range(len(opt_f))
    headers = ["opt_i", "opt_f", "opt_sol", "opt_comb_id",  "opt_intersect_id"]
    table = zip(opt_i, opt_f, opt_sol, opt_comb_id,  opt_intersect_id)
    print(tabulate(table, headers, tablefmt="fancy_grid"))
if p_i == 0:
    print("")
    print("")



####### Finding if the optimum points are feasible ############################################################################
###################################################################################################
###################################################################################################

# calculate or load from output of the previous run
if not load_feasible: # load feasible solutions from output of the previous run
    # This function checks if optimums found in the previous step are feasible. It first plugs the values of variables at the optimum
    # into constraints and then combines them with the constraints on the costant to make a system of equations that it solves for the
    # unknown constant and returns a feasible interval for the unknown constant. It can be empty for optimums that are always infeasible.
    #
    # Outputs:
    # feasible_interval: the interval in which the solution for optimum is feasible
    # feasible_f: the value of optimum f from the previous step
    # feasible_sol: values of variables at the optimum point from the previous step
    # feasible_comb_id:#  combination id from the previous step for backward check
    # feasible_intersect_id: intersection id from the previous step for backward check
    # feasible_opt_id: optimum id from the previous step for backward check


    feasible_interval, feasible_f, feasible_solution, feasible_combination = feasible(opt_f, opt_solution, opt_combination, Constants, c_exp, c2_exp, p_n, p_i)



    # free memory
    opt_f = None
    opt_solution = None
    opt_combination = None

    if p_n==1 and p_i==0:  # only save in non-parallel run (saving for troubleshooting)
        with open('./output/5_feasible.txt', 'wb') as fp:
            pickle.dump([feasible_interval, feasible_f, feasible_solution, feasible_combination], fp)
else:
    with open('./output/5_feasible.txt', 'rb') as fp:
        feasible_interval, feasible_f, feasible_solution, feasible_combination = pickle.load(fp)


# print all feasibles

if p_i == 0 and print_feasible:
    print("")
    print("List of feasible optimums:")
    feasible_i = range(len(feasible_f))
    headers = ["feasible_i", "feasible_f", "feasible_sol", "feasible_comb_id"  ]
    table = zip(feasible_i, feasible_f, feasible_solution, feasible_combination )
    print(tabulate(table, headers, tablefmt="fancy_grid"))
if p_i == 0:
    print("")
    print("")

# for i in range(len(feasible_solution)):
#     int = feasible_interval[i]
#     if str(int) != "FiniteSet(0)":
#         print(feasible_solution[i])
#         print(feasible_interval[i])
# exit()



####### Finding which feasible optimums are maximums ##############################################
###################################################################################################
###################################################################################################

# calculate or load maximums
if not load_max: # load maximums
    # in maximum_f function, we start with an array than contain one element, the full allowed interavl of the constant.
    # We also set the optimum function value for this iterval to be negative infinity (-sp.oo).
    # later we update the arrays to have smaller inetrvals each with the optimum value for that interval.
    # Outputs:
    # max_f: an array of maximum values for all the subintervals of possible interval of the constant
    # max_sol: to store the variable values for each interval
    # max_comb_id: keeping id for backward check
    # max_intersect_id: keeping id for backward check
    # max_opt_id: keeping id for backward check
    # max_feasible_id: keeping id for backward check
    # max_feasible_interval: list of intervals where the current maximum is feasible
    # (can be larger than interval below sincs below is both feasible and maximum)
    # max_interval: array of subintervals of the allowed region of the constant that has a distict maximum value save in max_f above.
    max_f, max_solution, max_combination, max_feasible_interval, max_interval = maximum_f(feasible_f, feasible_solution, feasible_combination, feasible_interval, c2_exp, comm, p_n, p_i)

    # in case of parallel processing, combining solutions from different processors
    max_f, max_solution, max_combination, max_feasible_interval, max_interval = collect_maximums(comm, p_n, p_i, max_f, max_solution, max_combination, max_feasible_interval, max_interval)

    if p_n==1 and p_i==0:  # only save in non-parallel run (saving for troubleshooting)
        with open('./output/6_maximum.txt', 'wb') as fp:
            pickle.dump([max_f, max_solution, max_combination, max_feasible_interval, max_interval], fp)
else:
    with open('./output/6_maximum.txt', 'rb') as fp:
        max_f, max_solution, max_combination, max_feasible_interval, max_interval = pickle.load(fp)




########################################################
############  sorting and printing results #############
########################################################


if p_i == 0:
    # finding middle of the intervals for sorting and saving to max_interval_mid
    max_interval_mid = [0] * len(max_interval)
    for i in range(len(max_interval)):
        i_s = list(max_interval[i].boundary)[0]
        i_end = list(max_interval[i].boundary)[-1]
        max_interval_mid[i] = (i_end + i_s )/2.0
        # evaluate functions for single point intervals
        if len(Constants) == 1 and i_end - i_s == 0:
            rep_dict = {str(Constants[0]):i_s}
            max_f[i] = max_f[i].subs(rep_dict)



    # sorting results based on middle of interval
    max_interval_mid, max_f, max_solution, max_combination,  max_feasible_interval, max_interval = (list(t) for t in zip(*sorted(zip(max_interval_mid, max_f, max_solution, max_combination,  max_feasible_interval, max_interval))))
    max_i = range(len(max_f))


    print(" ")
    print(" >>>> Final Results <<<<<< ")
    print(" ")
    print("Number of maximum intervals: ", len(max_f))
    print(" ")

    # printing combinations for each maximum in a table
    print("> List of combinations of constraints from which each maximum is found is saved to /output/0_max_out_2.txt")
    headers_2 = ["max_i",  "max_combination"]
    table_2 = zip(max_i,  max_combination)
    table_p_2 = tabulate(table_2, headers_2)

    # saving to file
    with open('./output/0_max_out_2.txt', 'w') as outputfile:
        outputfile.write("List of combinations of constraints from which the maximum is found:")
        outputfile.write("\n")
        outputfile.write(table_p_2)
        outputfile.write("\n")
        outputfile.write("max_i: id of the maximum. \n")
        outputfile.write("max_intersection_f: The function found as the result of intersection of combination of constraints and the main function. \n")
        outputfile.write("max_combination: combination of constraints from which the maximum is found. \n")




    # printing all maximums in a table
    print(" ")
    print("List of intervals of the unknown constant and the maximum values in those:")
    print("(also saved to /output/0_max_out_1.txt)")
    headers_1 = ["maximum id", "max interval", "max function value", "max input values"]
    table_1 = zip(max_i, max_interval, max_f, max_solution)
    table_p_1 = tabulate(table_1, headers_1)
    print(table_p_1)
    print("-----------------------------------")
    print("function:", f)
    print("-----------------------------------")
    print("maximum id: id of the solution in the output list of maximum function()")
    print("max interval: The intervals of the unknown constant over which the answer is valid.")
    print("(Just shows {0} if there is no unknown constant.)")
    print("max function value: The function value of the maximum (can be dependent on the unknown constant).")
    print("max input values: The values of the input variables at the maximum (can be dependent on the unknown constant).")
    # saving to file
    with open('./output/0_max_out_1.txt', 'w') as outputfile:
        outputfile.write("List of intervals of the unknown constant and the maximum values in those: \n")
        outputfile.write(table_p_1)
        outputfile.write("\n")
        outputfile.write("-----------------------------------\n")
        outputfile.write("function: "+ str(f)+"\n")
        outputfile.write("-----------------------------------\n")
        outputfile.write("\n")
        outputfile.write("maximum id: id of the solution in the output list of maximum function() \n")
        outputfile.write("max interval: The intervals of the unknown constant over which the answer is valid. \n")
        outputfile.write("(Just shows {0} if there is no unknown constant.) \n")
        outputfile.write("max function value: The function value of the maximum (can be dependent on the unknown constant). \n")
        outputfile.write("max input values: The values of the input variables at the maximum (can be dependent on the unknown constant). \n")











