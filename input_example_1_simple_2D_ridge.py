### To run the code for this example copy contents of this
# example to "input.py" and run "main.py".

### simple 2D quadratic function with a ridge

####### Input Function and Constraints ############################################################################
############################################################################################

# main function. A string in latex format.
Function = '2-(x_1+y_1)^2'

# constraints on independent variables.
# Add a constraint to the list using Constraints.append() command.
# The constraint needs to be a string in latex format.
# Right hand side of inequality must always be zero.
Constraints = []
Constraints.append('x_1 \geq 0' )
Constraints.append('x_1 - 1 \leq 0')
Constraints.append('y_1 \geq 0')
Constraints.append( 'y_1 - 1 \leq 0')
Constraints.append('x_1 + y_1 - s \geq 0' )

# list of constants. The constant name can only be o, s or m.
# Only one constant is supported for this version.
Constants = [s]

# constraints on constants
Constants_Constraints = []
Constants_Constraints.append('s \geq 0')
Constants_Constraints.append('s - 2 \leq 0')

