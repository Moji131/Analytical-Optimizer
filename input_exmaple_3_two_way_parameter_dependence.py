### To run the code for this example copy contents of this
# example to "input.py" and run "main.py".

# function and constraints for CHSH inequality with two-way parameter dependence.

####### Input Function and Constraints ############################################################################
############################################################################################

# main function. A string in latex format.
Function = ' (1 + 4m_1n_1 - 2m_1 - 2n_1) + (1 + 4m_2n_2 - 2m_2 - 2n_2) + (1 + 4m_3n_3 - 2m_3 - 2n_3) - (1 + 4m_4n_4 - 2m_4 - 2n_4) '

# constraints on independent variables.
# Add a constraint to the list using Constraints.append() command.
# The constraint needs to be a string in latex format.
# Right hand side of inequality must always be zero.

Constraints = []
Constraints.append('m_1-m_2 - s \leq 0' )
Constraints.append('m_2-m_1 - s \leq 0 ')
Constraints.append('m_3-m_4 - s \leq 0 ')
Constraints.append('m_4-m_3 - s \leq 0 ')

Constraints.append('n_1-n_3 - s \leq 0 ')
Constraints.append('n_3-n_1 - s \leq 0 ')
Constraints.append('n_2-n_4 - s \leq 0 ')
Constraints.append('n_4-n_2 - s \leq 0 ')

Constraints.append('m_1 \geq 0')
Constraints.append('m_1 - 1 \leq 0')
Constraints.append('m_2 \geq 0')
Constraints.append('m_2 - 1 \leq 0')
Constraints.append('m_3 \geq 0')
Constraints.append(r'm_3 - 1 \leq 0')
Constraints.append('m_4 \geq 0')
Constraints.append('m_4 - 1 \leq 0')

Constraints.append('n_1 \geq 0')
Constraints.append('n_1 - 1 \leq 0')
Constraints.append('n_2 \geq 0')
Constraints.append('n_2 - 1 \leq 0')
Constraints.append('n_3 \geq 0')
Constraints.append(r'n_3 - 1 \leq 0')
Constraints.append('n_4 \geq 0')
Constraints.append('n_4 - 1 \leq 0')

# list of constants. The constant name can only be o, s or m.
# Only one constant is supported for this version.
Constants = [s]

# constraints on constants
Constants_Constraints = []
Constants_Constraints.append('s \geq 0')
Constants_Constraints.append('s - 1 \leq 0')

