####### options ############################################################################
############################################################################################

# n single processor running, at the and of each step the results can be saved to a file.
# By setting load_<step> = False the code calculates normally and
# overwrites the previously saved file.
# If set to True it skips the calculation and loads the previous results from file
# This is useful if one wants to print out and check outputs of each stage without redoing the calculations.

load_input =  False
load_constraint_combinations =  False
load_intersect =  False
load_optimum =  False
load_feasible =  False
load_max = False

# load_input = True   # Always load input if you are loading any of the below. Translation between working variables and input variables chagnes everytime that you read from input.py
# load_constraint_combinations = True
# load_intersect = True
# load_optimum = True
# load_feasible = True
# load_max = True

### Whether or not to print results of initial stages.
# These can be very long prints for large functions and high number of constraints
print_intersects = False
print_optimum = False
print_feasible = False

# print_intersects = True
# print_optimum = True
# print_feasible = True




