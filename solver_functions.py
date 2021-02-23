import sympy as sp
from sympy.solvers import solve, nonlinsolve, solveset, linsolve
import time
import datetime
import MPI_functions as mpi



####### Defining optimization functions ############################################################################
###################################################################################################
###################################################################################################



def intersect(f, constraint_combinations, var_list, p_i):
    # This function calculates the intersection of the function f with each combination
    # of constraints in the list "comb".

    # intersect_f: list to save intersection of f with set of constraints
    # intersect_sol: list to save intersection of constraints with each other
    # intersect_comb_id: to keep id of combination for backward check


    intersect_f=[] # to keep intersction function
    intersect_solution= [] # to keep solution to intersection of combination of constraints
    intersect_combination = [] # to kkep combination id

    # saving time for timing the loop
    # start = time.time()
    # t_p = t_print # setting first progress printing time


    for comb_id in range(len(constraint_combinations)):
        #this loop gets the constraint combinations (inequlity must be turned into equality),
        # and finds theire intersection with the main function

        #########################################################################################################
        ### main part of the caulation is to use solve from Sympy to solve the sytem of equations in comb[comb_id]
        # and find values for variables in Q_list
        sol = solve(constraint_combinations[comb_id], var_list)  # finds intersection of constraints set
        sol_n = len(sol) # saving the number of solutions

        if (sol_n != 0):   #if intersection is not empty add the solution to array
            # replacing indexed variables in f with their values in the solution to evaluate function value.
            # Might still not be a number and depend on some of the variables
            intersect_f = intersect_f + [f.subs(sol)]
            intersect_solution = intersect_solution + [sol]   # adding solution of constraint set to an array
            intersect_combination = intersect_combination + [constraint_combinations[comb_id]]  # keeping id of constraint combination for backward check

        # printing the progress by processor 0
    #
    #     t = time.time() - start
    #     if t > t_p :
    #         perc = ((comb_id+1)/(len(constraint_combinations)-1)*100)
    #         t_left = t * (100 - perc) / perc
    #         print("process",p_i, "intersect", int(perc), "%  Time to finish:" , str(datetime.timedelta(seconds=int(t_left)) ) )
    #         t_p = t_p + t_print
    #
    # # calculating and printing total time to calculate
    # end = time.time()
    # print("process",p_i," Intersect completed. Total time: ",  str(datetime.timedelta(seconds=(end-start))) )
    print("process", p_i, " Intersect completed.")


    # intersect_f: list to save intersection of f with set of constraints
    # intersect_sol: list to save intersection of constraints with each other
    # intersect_comb_id: to keep id of combination for backward check
    return intersect_f, intersect_solution, intersect_combination



#################################
#################################
#################################



def optimum(intersect_f, intersect_solution, intersect_combination, var_list, p_n, p_i):
    #This function takes the partial derivativs of the all intersections of function f with combinations of constraints
    # with respect to all the variables usinf diff() function of Sympy then it solves the system of equations of derivatives
    # for zeros. It outputs a solution list (values of variables at optimum): "sol_list"
    # and corresponding optimum function values list: "f_opt".

    # Outputs:
    # opt_f: optimum value of f intersection with constraint combination (optimum values withing the volume of and at the boundaries)
    # opt_sol: to keep values of variables for each optimum
    # opt_comb_id: to keep combination id for backward check
    # opt_intersect_id: to keep intersection id for backward check

    var_n = len(var_list)

    opt_f = [] # to keep optimum value of f (intersection with constraint set)
    opt_solution = [] # to keep values of variables for this optimum (from previous step)
    opt_combination = [] # to keep combination id (from previous step)

    # start = time.time() # used for timing
    # t_p = t_print # setting the first time to print


    # looping over own portion of intersection arrays
    for i in range(len(intersect_f)):
        #making a list of partial derivatives with respect to all variables
        df=[0 for i in range(var_n)] #empty list to store partial derivatives
        for var_i in range(var_n):   #loop over all vraibles
            df[var_i]=sp.diff(intersect_f[i],var_list[var_i])  #caculating partial derivative


        test_df = [0 for i in range(var_n)]
        if df == test_df : # if function is constant (all derivates are zero) just add the function to list
            if len(intersect_solution[i]) == var_n: # if the solution is fully determined means intersection is a point (if it is not fully detemined it has a ridge. we do not save ridge and it appears as a point in one of the other boundaries )
                opt_f = opt_f + [intersect_f[i]]
                opt_solution = opt_solution + [intersect_solution[i]]
                opt_combination = opt_combination + [intersect_combination[i]]

        else: # if f is not constant
            # solving system of equations of derivatives for zeros
            ###### This is the heart of optimization. The main limitaion of the code is the ability of nonlinsolve() of Sympy
            # to sovle system of nonlinear equations. If the main function is quadratic this sytem is linear and there is certainly
            # no problem. For nonlinear functions more complicated than quadratic it may get stuk here and run forever.
            sol = nonlinsolve(df, var_list)
            sol_n = len(sol) # saving the number of sulutions
            if sol_n != 0 :
                for sol_i in range(sol_n):  # iteration over all the zero derivative solutions
                    # making the optimal solution by combining the solution of intersection (values of variables)
                    # and optimization solutions
                    opt_sol_dict = dict(zip(var_list, list(sol)[sol_i]))
                    f1_list = []
                    for k, v in opt_sol_dict.items():
                        f1_list = f1_list + [k-v]
                    for k, v in intersect_solution[i].items():
                        f1_list = f1_list + [k-v]
                    sol_all = nonlinsolve(f1_list,var_list)
                    dict_all = dict(zip(var_list, list(sol_all)[0]))
                    opt_sol_o = dict_all

                    # check if we have ridges (solution is underdetermined)
                    no_ridge = True
                    for k, v in opt_sol_o.items():
                        if k == v: # check if the solution for a variable is not a number but itself (means it is a free variable)
                            no_ridge = False

                    if no_ridge: # throwing away ridges
                        opt_f_o = intersect_f[i]
                        opt_f_o = sp.simplify(opt_f_o.subs(opt_sol_o))  #replacing indexed variables with their values in the solution to evaluate function value

                        opt_f = opt_f + [opt_f_o] # to keep optimum value of f (intersection with constraint set)
                        opt_solution = opt_solution + [opt_sol_o] # to keep variable values at the optimum points of f (intersection with constraint set)
                        opt_combination = opt_combination + [intersect_combination[i]] # to keep combination id (for backward check)

    #     # print timing and progress
    #     t = time.time() - start
    #     if t > t_p :
    #         perc = ((i+1)/(len(intersect_f)-1)*100)
    #         t_left = t * (100 - perc) / perc
    #         print("process",p_i, "optimum", int(perc), "%  Time to finish:" , str(datetime.timedelta(seconds=int(t_left)) ) )
    #         t_p = t_p + t_print
    #
    # # calculate and print total time for optimum calculation
    # end = time.time()
    # print("process", p_i, " Optimum completed. Total time: ",  str(datetime.timedelta(seconds=(end-start))) )

    print("process", p_i, " Optimum completed.")

    # Outputs:
    # opt_f: optimum value of f intersection with constraint combination (optimum values withing the volume of and at the boundaries)
    # opt_sol: to keep values of variables for each optimum
    # opt_comb_id: to keep combination id for backward check
    # opt_intersect_id: to keep intersection id for backward check

    return opt_f, opt_solution, opt_combination



#################################
#################################
#################################


def feasible(opt_f, opt_solution, opt_combination, Constants, c_exp, c2_exp, p_n, p_i):
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
    # feasible_opt_id: optimum id from the previous step for backward



    feasible_interval = []
    feasible_f = []
    feasible_solution = []
    feasible_combination = []

    # start = time.time() # is used to calculate total time
    # t_p = t_print # setting the first time to print progress

    # going over all optimum points
    for i in range(len(opt_f)):

        # making a system of equations of constarints on variables and constants
        sys_ineq = []
        for j in range(len(c_exp)):
            sys_ineq = sys_ineq + [ c_exp[j].subs(opt_solution[i]) ] # substituting the optimum point in constarints on variables and adding to system
        for j in range(len(c2_exp)):
            sys_ineq = sys_ineq + [ c2_exp[j] ] # adding constraints on constants to the system


        sol = solve(sys_ineq, Constants, set=True) # solving the system of equation for unknown constants

        dic1 = opt_solution[i]
        for key, val in dic1.items():
            if str(val) == "-s" and not (sol == False) and not ( str(sol) == "Eq(s, 0)"):
                print(dic1)
                print("sys_ineq", sys_ineq)
                print("sol", sol)
                input("enter")


        # print(p_i, constants)

        if sol != False: # throw away optimum points that are always infeasible

            if sol == []: # for the case that there is no constant just save a dummy interval
                sol_set = sp.Interval(0,0)
            else: # saving interval of feasibility as a set
                sol_set = sol.as_set()




            # appending to ouput arrays
            feasible_interval = feasible_interval + [sol_set] # inetrval of feasibility
            feasible_f = feasible_f + [opt_f[i]] # optimum value of f
            feasible_solution = feasible_solution + [opt_solution[i]] # values of variables for the optimum value
            feasible_combination = feasible_combination + [opt_combination[i]]



        # printing progress
    #     t = time.time() - start
    #     if t > t_p :
    #         perc = ((i+1)/(len(opt_f)-1)*100)
    #         t_left = t * (100 - perc) / perc
    #         print("process",p_i, "feasible", int(perc), "%  Time to finish:" , str(datetime.timedelta(seconds=int(t_left)) ) )
    #         t_p = t_p + t_print
    #
    #
    # # calculating and printing total time
    # end = time.time()
    # print("process", p_i, " Feasible completed. Total time: ",  str(datetime.timedelta(seconds=(end-start))) )
    print("process", p_i, " Feasible completed.")


    # Outputs:
    # feasible_interval: the interval in which the solution for optimum is feasible
    # feasible_f: the value of optimum f from the previous step
    # feasible_sol: values of variables at the optimum point from the previous step
    # feasible_comb_id:#  combination id from the previous step for backward check
    # feasible_intersect_id: intersection id from the previous step for backward check
    # feasible_opt_id: optimum id from the previous step for backward check
    return feasible_interval, feasible_f, feasible_solution, feasible_combination



#################################
#################################
#################################



# finds maximum for different intervals given feasible answers
def maximum_f(feasible_f, feasible_solution, feasible_combination, feasible_interval, c2_exp, comm, p_n, p_i):
    # this function gets the array of optimum values and their feasible regions and comparing all of them finds
    # the maximum value for diffrent inetrvals of possible reagion of the constant.
    # the maximum value can be a number or a function of the constant.


    if c2_exp == []: # if there is no constants
        sol1_set = sp.Interval(0,0) # set inetval
    else:
        # convert constarint on constants to set. This is the initial interval to be divided later for diffrent maximums
        sol1 = sp.solve(c2_exp)
        sol1_set = sol1. as_set()

    # here we define an array of intervals that have one interavl which is allowed region of the constant
    # We also set the optimum function value for this iterval to be negative infinity (-sp.oo)
    # later we update the arrays to have smaller inetrvals each with the optimum value for that interval
    # Outputs:
    # max_f: an array of maximum values for all the subintervals of possible interval of the constant
    # max_solution: to store the variable values for each interval
    # max_comb_id: keeping id for backward check
    # max_intersect_id: keeping id for backward check
    # max_opt_id: keeping id for backward check
    # max_feasible_id: keeping id for backward check
    # max_feasible_interval: list of intervals where the current maximum is feasible
    # (can be larger than interval below sincs below is both feasible and maximum)
    # max_interval: array of subintervals of the allowed region of the constant that has a distict maximum value save in max_f above.



    max_f = [-sp.oo] # set the maximum value of function to negative infinity
    max_solution = ["No answer"] # to store the variable values for each maximum
    max_combination = ["No answer"]
    max_feasible_interval = ["No answer"]
    max_interval = [sol1_set] # set array of interval to have one member which is the allowed reagion of constant


    # t_p = t_print
    # # used to caculate total caculate time
    # start = time.time()

    # loop over all feasible solutions
    for i_f in range(len(feasible_f)):
        # max_2* arrays are used as temporary storages of maximum intervals and at the end of the loop are
        # copied to max* arrays
        max2_f = []
        max2_sol = []
        max2_combination = []
        max2_feasible_interval = []
        max2_interval = []

        # loop over all previous intervals of maximum values saved in max* arrays
        for i_max in range(len(max_interval)):
            # intersection of current optimum interval with each previously calculated maximum interval
            inter1 = sp.Intersection(max_interval[i_max], feasible_interval[i_f])

            if inter1 == sp.EmptySet: # if there is no intersection retain the current interval and associated maximum
                max2_f.append(max_f[i_max])
                max2_sol.append(max_solution[i_max])
                max2_combination.append(max_combination[i_max])
                max2_feasible_interval.append(max_feasible_interval[i_max])
                max2_interval.append(max_interval[i_max])
            else: # if there is an intersection, find the interval that the new maximum is larger than the previous maximum
                exp1 = (feasible_f[i_f] - max_f[i_max]) > 0 # make an expersion to solve to see if current optimum is larger than previously found optimum for this interval
                exp1 = exp1.subs({'s': 'd'})  # becuase sympy is stupid and sometimes considers same 's' two variables.
                exp1 = exp1.subs({'d': 's'})
                if exp1 == False: # new optimum is not larger anywhere
                    sol1_set = sp.EmptySet
                elif exp1 == True: # new interval is always bigger
                    sol1_set = sp.UniversalSet
                else: # if only larger on some interval, find the interval
                    sol1 = solve(exp1)
                    sol1_set = sol1.as_set()

                new_max_inter = sp.Intersection(inter1, sol1_set) # define a new interval that is intersection of where the new optimum is larger and is feasible and whithin current maximum interval
                if new_max_inter != sp.EmptySet: # if this new interval is not empty append this interval and its f value to the list of intervals for maximum
                    max2_f.append(feasible_f[i_f])
                    max2_sol.append(feasible_solution[i_f])
                    max2_combination.append(feasible_combination[i_f])
                    max2_feasible_interval.append(feasible_interval[i_f])
                    max2_interval.append(new_max_inter)

                    max_comp_inter = sp.Complement(max_interval[i_max], new_max_inter) # find the complement of the inetvarl above to retain the previous maxium value for the remaining interval
                    if max_comp_inter != sp.EmptySet: # if remaining is not empty set add it to the array of maximums and retain the previous value for the maximum
                        max2_f.append(max_f[i_max])
                        max2_sol.append(max_solution[i_max])
                        max2_combination.append(max_combination[i_max])
                        max2_feasible_interval.append(max_feasible_interval[i_max])
                        max2_interval.append(max_comp_inter)

                else: # if this new interval is empty retain the previous interval and its value as maximum
                    max2_f.append(max_f[i_max])
                    max2_sol.append(max_solution[i_max])
                    max2_combination.append(max_combination[i_max])
                    max2_feasible_interval.append(max_feasible_interval[i_max])
                    max2_interval.append(max_interval[i_max])

        # copy from temporary to permanent storage for the next iteration
        max_f = max2_f.copy()
        max_solution = max2_sol.copy()
        max_combination = max2_combination.copy()
        max_feasible_interval = max2_feasible_interval.copy()
        max_interval = max2_interval.copy()

        # print progress

    #     t = time.time() - start
    #
    #     # printing progress
    #     if t > t_p :
    #         perc = ((i_f+1)/(len(feasible_f)-1)*100)
    #         if perc > 0:
    #             t_left = t * (100 - perc) / perc
    #             print("process",p_i, "maximum", int(perc), "%  Time to finish:" , str(datetime.timedelta(seconds=int(t_left)) ) )
    #         t_p = t_p + t_print
    #
    # # caculate and print run time
    # end = time.time()
    #
    # print("process", p_i, " Maximum completed. Total time: ",  str(datetime.timedelta(seconds=(end-start))) )
    print("process", p_i, " Maximum completed.")
    # Outputs:
    # max_f: an array of maximum values for all the subintervals of possible interval of the constant
    # max_solution: to store the variable values for each interval
    # max_combination: keeping initial combination of constraints that lead to this.
    # max_feasible_interval: list of intervals where the current maximum is feasible
    # (can be larger than interval below sincs below is both feasible and maximum)
    # max_interval: array of subintervals of the allowed region of the constant that has a distict maximum value save in max_f above.
    return max_f, max_solution, max_combination, max_feasible_interval, max_interval


#################################
#################################
#################################




def maximum_2(max_f, max_solution, max_combination, max_feasible_interval, max_interval, max3_f, max3_solution, max3_combination, max3_feasible_interval, max3_interval):
    # this function gets two arrays of maximum intervals and returns an overall maximum interval array
    # the concept is simillar to the function maximum_f() above. Only the starting point is the first
    # array, max_interval (instead the total feasible reagion of the constant ),
    # and the loop is over the second array (max3_interval)
    # this is used to combine intervals and maximums caculated by diffrent processors in case of parallel processing.

    # loop over all inetrvals in max3_interval
    for i_max3 in range(len(max3_interval)):

        # temporary storage for intervals and their values
        max2_f = []
        max2_solution = []
        max2_combination = []
        max2_feasible_interval = []
        max2_interval = []

        # loop over all inetrvals in max_interval to be updated by max3_interval
        for i_max in range(len(max_interval)):
            # intersection of current maximum interval with each new maximum interval
            inter1 = sp.Intersection(max_interval[i_max], max3_interval[i_max3])
            if inter1 == sp.EmptySet: # if there is no intersection keep current value
                max2_f.append(max_f[i_max])
                max2_solution.append(max_solution[i_max])
                max2_combination.append(max_combination[i_max])
                max2_feasible_interval.append(max_feasible_interval[i_max])
                max2_interval.append(max_interval[i_max])
            else: # if there is intesection
                # check in which interval the new maximum is larger than the current one
                exp1 = sp.sympify(max3_f[i_max3] - max_f[i_max] > 0)
                exp1 = exp1.subs({'s': 'd'})  # becuase sympy is stupid and sometimes considers same 's' two variables.
                exp1 = exp1.subs({'d': 's'})
                if exp1 == False:
                    sol1_set = sp.EmptySet
                elif exp1 == True:
                    sol1_set = sp.UniversalSet
                else:
                    sol1 = solve(exp1)
                    sol1_set = sol1.as_set()

                new_max_inter = sp.Intersection(inter1, sol1_set)  # intersection of interval that the new one is bigger and feasible
                if new_max_inter != sp.EmptySet: # if not empty make a new interval and keep the new maximum value as its maximum
                    max2_f.append(max3_f[i_max3])
                    max2_solution.append(max3_solution[i_max3])
                    max2_combination.append(max3_combination[i_max3])
                    max2_feasible_interval.append(max3_feasible_interval[i_max3])
                    max2_interval.append(new_max_inter)

                    max_comp_inter = sp.Complement(max_interval[i_max], new_max_inter) # complement of the updated interval
                    if max_comp_inter != sp.EmptySet: # if complemet is not emty, ratain the previous maximum value for this complement interval
                        max2_f.append(max_f[i_max])
                        max2_solution.append(max_solution[i_max])
                        max2_combination.append(max_combination[i_max])
                        max2_feasible_interval.append(max_feasible_interval[i_max])
                        max2_interval.append(max_comp_inter)
                else: # if there is no intersection between current interval and the interval that the new one is bigger ratin the current value
                    max2_f.append(max_f[i_max])
                    max2_solution.append(max_solution[i_max])
                    max2_combination.append(max_combination[i_max])
                    max2_feasible_interval.append(max_feasible_interval[i_max])
                    max2_interval.append(max_interval[i_max])

        # move from temporary storage to permanent
        max_f = max2_f.copy()
        max_solution = max2_solution.copy()
        max_combination = max2_combination.copy()
        max_feasible_interval = max2_feasible_interval.copy()
        max_interval = max2_interval.copy()

    return max_f, max_solution, max_combination, max_feasible_interval, max_interval




