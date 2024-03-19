#Computational Genomics Spring 2023 Final Project
#Jamie Bregman - 5/10/2023

"""
Python script that optimizes 15 parameters for the calculation of the p-ratio, or the retention of genes after consecutive whole genome duplication events.
This software uses "grid" optimization, meaning that each value for each parameter is tested against one another to determine the optimized values.
Requires known values for time of the first WGD event, time of the second WGD event, and the known p-ratio for those events.
Optimizes the parameters based on equations from Amanda and Konrad et Al.


3 types of gene outcome that were optimized based on percentages of the genome, meaning that these 3 parameters add up to 1:
    alt-functionalization: genes can sub or neo functionalize
        - alpha_alt
    dosage balance: genes are sensitive to dosage balance effects
        - alpha_dos
    non-functionalization: both gene copies can only be retained by chance
        - alpha_non
    
4 parameters were optimized for EACH type of gene outcome (sub/alt-functionalization, dosage balance, non-functionalization).
    b,c: shape of the "curve", i.e., how the dynamics behave when moving from instaneous rate to asymptotic rate.
        - alt_b, alt_c, 
        - dos_b, dos_c
        - non_b, non_c
    d,f: rate at which fully redundant genes get lost from the genome (d+f diminishes to d, which is the rate at which non-duplicated genes are lost from the genome).
        - alt_d, alt_f
        - dos_d, dos_f
        - non_d, non_f

This software can be helpful when modelling gene retention for given WGD events and observing how parameters change depending on the times of these events.
Additionally, because the parameters given in Konrad et Al. have such wide bounds, this software allows you to better understand the way(s) in which the equation behaves for a given set of WGD events.
"""

#importing modules
import math
import numpy as np
import itertools
from scipy.optimize import minimize


time1 = float(input("Please enter time of first WGD event (time 1): "))
time2 = float(input("Please enter time of second WGD event (time 2): "))
known_p_ratio = float(input("Please enter the p-ratio of gene retention for these WGD events: "))

"""
#Known values example
time1 = 0.49
time2 = 0.75
known_p_ratio = 1.08
"""

#Survival function describing the corresponding probability of survival to a time used in calculating p-ratio (taken from Amanda's code / Konrad et. Al)
#Function uses "try" and "except" as certain parameter combinations result in math Overflow Errors.
#However, the parameter bounds could not be changed to fix these errors without getting rid of parameter values that may work with other combinations of parameters.
#Thus, the try and except method was employed.

def prob_surv(b, c, d, f, time):
    n_max = 100 #n_max value taken from Konrad et. Al
    summation = 0
    for n in range(0, n_max):
        nfac = math.factorial(n)
        beta = (((-b) ** n) * (time ** ((c * n) + 1))) / ((c * n * nfac) + nfac)
        summation += beta
    try:
        survival_probability = math.exp(-d * time - f * summation)
    except OverflowError:
        return None
    return survival_probability if not math.isnan(survival_probability) else None


#Function for calculating the p-ratio given the proper parameters (taken from Amanda's paper / Konrad et. Al)
#This function uses if statements to return None if any of the parameters result in an OverflowError when calculating the function above.
#Additionally, this function only returns values when alpha_alt and alpha_dos are less than or equal to 1
    #This ensures that the genome percentages are viable values. Otherwise, due to the range of parameter combinations and the "grid" aspect, this function would yield values that have non-viable values for the alpha values.
def pratio_function(b_alt, c_alt, d_alt, f_alt, b_dos, c_dos, d_dos, d_non, f_non, alpha_alt, alpha_dos, time1, time2):

    if alpha_alt + alpha_dos <= 1:
    
        f_dos = d_dos * -1
        alpha_non = 1 - alpha_alt - alpha_dos

        s_alt_time1 = prob_surv(b_alt,c_alt,d_alt,f_alt,time1)
        if s_alt_time1 is None:
            return None
        
        s_alt_time2 = prob_surv(b_alt,c_alt,d_alt,f_alt,time2)
        if s_alt_time2 is None:
            return None
            
        s_dos_time1 = prob_surv(b_dos,c_dos,d_dos,f_dos,time1)
        if s_dos_time1 is None:
            return None
        
        s_dos_time2 = prob_surv(b_dos,c_dos,d_dos,f_dos,time2)
        if s_dos_time2 is None:
            return None

        #b_non and c_non are always 0 and 1, respectively
        s_non_time1 = prob_surv(0,1,d_non,f_non,time1)
        if s_non_time1 is None:
            return None
        
        s_non_time2 = prob_surv(0,1,d_non,f_non,time2)
        if s_non_time1 is None:
            return None
        
        num = (2*alpha_alt*s_alt_time1*s_alt_time2) + (2*alpha_dos*s_dos_time1*s_dos_time2) + (2*alpha_non*s_non_time1*s_non_time2)
        denom = ((1-s_alt_time1)*alpha_alt*s_alt_time2) + ((1-s_dos_time1)*alpha_dos*s_dos_time2) + ((1-s_non_time1)*alpha_non*s_non_time2)
        temp = num/denom

        num1 = ((1-s_alt_time1)*alpha_alt) + ((1-s_dos_time1)*alpha_dos) + ((1-s_non_time1)*alpha_non)
        denom1 = (2*alpha_alt*s_alt_time1) + (2*alpha_dos*s_dos_time1) + (2*alpha_non*s_non_time1)
        temp1 = num1/denom1

        pratio = temp * temp1
        
        return pratio

    else:
        return None


#Defining the bounds of the paramaters
#Bounds were created using Konrad et. Al, refined using the example values in Table 2 from Konrad et. Al, and discussed with professor Liberles.
#'np.linspace' is used to create a NumPy array of evenly spaced values over a given interval.
    #The first 2 values in the function call indicate the interval in which to build the array over.
    #The last value in the function call indicates the number of evenly spaced intervals to create.
        #This value can be increased to increase the number of intervals and thus the precision of the optimization.
        #However, my machine could only handle 3 intervals per parameter without taking many hours to run.

#Parameters for neo-functionalization / sub-functionalization
b_alt_range = np.linspace(0.0001,100,3)  
c_alt_range = np.linspace(1,3,3)
d_alt_range = np.linspace(0.0001,20,3)
f_alt_range = np.linspace(0.0001,20,3)

#Parameters for Dosage
b_dos_range = np.linspace(-100,-0.0001,3)
c_dos_range = np.linspace(0.0001,0.9999,3)
d_dos_range =  np.linspace(-0.001,-0.000001,3)

#Parameters for Non-functionalization       
d_non_range = np.linspace(10,20,3)  
f_non_range = np.linspace(0.0001,20,3)

#Paramaters for the proportion of the initial genome for each category above
#These values go up to 20 intervals to ensure more percentage combos
alpha_alt_range = np.linspace(0.01,0.99,20)
alpha_dos_range = 1 - alpha_alt_range


#Grid search optimization
#Creating a list of arrays for each parameter bound
param_ranges = [b_alt_range, c_alt_range, d_alt_range, f_alt_range, b_dos_range, c_dos_range, d_dos_range, d_non_range, f_non_range, alpha_alt_range, alpha_dos_range]

#Using itertools to create each combination of each parameter (grid)
param_combinations = list(itertools.product(*param_ranges)) 

#Initializing the starting "minimum error"
lowest_error = float("inf")

#Initializing the starting "best parameters"
best_params = None

#Looping over each parameter combination and calculating the p-ratio
for params in param_combinations:
    temp_pratio = pratio_function(*params, time1, time2)

    #If that p-ratio is not None (i.e., if the parameter values do not cause an OverflowError)
    #Calculate the error of that p-ratio compared to the known value
    #Assign that new p-ratio to have the lowest_error and store the parameters for that value
    #Continue until the lowest error is found
    if temp_pratio is not None:
        error = abs(temp_pratio - known_p_ratio)
        if error < lowest_error:
            lowest_error = error
            best_params = params
    else:
        continue

"""
##alternatively, you can use nested for loops to create "grid" parameter combinations (less efficient)

def grid_optimize(b_alt_range, c_alt_range, d_alt_range, f_alt_range, b_dos_range, c_dos_range, d_dos_range, d_non_range, f_non_range, alpha_alt_range, alpha_dos_range, time1, time2, known_p_ratio):

    best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos = None, None, None, None, None, None, None, None, None, None, None #Initializing each optimized parameter as 'None'

    #Initializing the minimum error as infinity
    min_error = float("inf") 

    #Nested for loops that act as a "grid"
    #This means that each value for each parameter is tested against each other to find the optimal values
    for b_alt in b_alt_range:
        for c_alt in c_alt_range:
            for d_alt in d_alt_range:
                for f_alt in f_alt_range:
                    for b_dos in b_dos_range:
                        for c_dos in c_dos_range:
                            for d_dos in d_dos_range:
                                for d_non in d_non_range:
                                    for f_non in f_non_range:
                                        for alpha_alt in alpha_alt_range:
                                            for alpha_dos in alpha_dos_range:

                                                #Calculate the p-ratio for a given set of parameters 
                                                temp_pratio = pratio_function(b_alt,c_alt,d_alt,f_alt,b_dos,c_dos,d_dos,d_non,f_non,alpha_alt,alpha_dos,time1,time2)

                                                if temp_pratio is not None:

                                                    #The "error" of that p-ratio is its absolute difference compared to the known p-ratio
                                                    error = abs(temp_pratio - known_p_ratio)

                                                    #If that error is less than the minimum error, then that p-ratio (and its associated parameters) becomes the new minimum error
                                                    if error < min_error:
                                                        min_error = error
                                                        best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos = b_alt,c_alt,d_alt,f_alt,b_dos,c_dos,d_dos,d_non,f_non,alpha_alt,alpha_dos

                                                else:
                                                    continue
                                                
    #Returning the parameter values that return the value closest to the known p-ratio                                        
    return best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos

#Calling the function and assigning each optimized parameter value to an output variable
best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos = grid_optimize(b_alt_range, c_alt_range, d_alt_range, f_alt_range, b_dos_range, c_dos_range, d_dos_range, d_non_range, f_non_range, alpha_alt_range, alpha_dos_range, time1, time2, known_p_ratio)
"""


#Assigning the optimized parameter values
best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos = best_params
calc_p_ratio = pratio_function(best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos, time1, time2)

#Creating a final list of optimized parameter values and adding in the parameter values that were not "bounded" (i.e., did not need to be 'optimized')
optimized_parameter_list = [best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos]
optimized_parameter_list.insert(7, optimized_parameter_list[6] * -1) #dosage balance "f" is equal to negative dosage balance "d"
optimized_parameter_list.insert(8, 0.0) #non-functionalization "b" is always equal to 0.0
optimized_parameter_list.insert(9, 1.0) #non-functionalization "c" is always equal to 1.0
optimized_parameter_list.append(1 - optimized_parameter_list[12] - optimized_parameter_list[13]) #non-functionalization genome percentage is 1 - altfunctionalization percentage - dosage balance percentage)

#Writing results to a file
with open("optimization_output.txt", "a") as file:

    #Writing the known values for time 1, time2, and p-ratio
    file.write("Optimized parameters for time 1 (" + str(time1) + ") time 2 (" + str(time2) + ") and p-ratio (" + str(known_p_ratio) + "):" + "\n")
    
    #Writing each parameter and its optimized value
    file.write("b_alt = " + str(optimized_parameter_list[0]) + "\n")
    file.write("c_alt = " + str(optimized_parameter_list[1]) + "\n")
    file.write("d_alt = " + str(optimized_parameter_list[2]) + "\n")
    file.write("f_alt = " + str(optimized_parameter_list[3]) + "\n")
    file.write("b_dos = " + str(optimized_parameter_list[4]) + "\n")
    file.write("c_dos = " + str(optimized_parameter_list[5]) + "\n")
    file.write("d_dos = " + str(optimized_parameter_list[6]) + "\n")
    file.write("f_dos = " + str(optimized_parameter_list[7]) + "\n")
    file.write("b_non = " + str(optimized_parameter_list[8]) + "\n")
    file.write("c_non = " + str(optimized_parameter_list[9]) + "\n")
    file.write("d_non = " + str(optimized_parameter_list[10]) + "\n")
    file.write("f_non = " + str(optimized_parameter_list[11]) + "\n")
    file.write("alpha_alt = " + str(optimized_parameter_list[12]) + "\n")
    file.write("alpha_dos = " + str(optimized_parameter_list[13]) + "\n")
    file.write("alpha_non = " + str(optimized_parameter_list[14]) + "\n")

    #Writing the resulting p-ratio achieved from these optimized values
    file.write("\n" + "Resulting p-ratio: " + str(calc_p_ratio) + "\n")

    #Writing the percent error between the "optimized" p-value and the known p-value
    file.write("Percent error: " + str(((abs(calc_p_ratio - known_p_ratio))/known_p_ratio)*100))


########################################################################################################################################################################################################################################################################################


#Fine-tuning the optimization off grid
#Using scipy minimize to fine-tune the hyper parameters to get more exact values
#This results in a smaller percent error and a value closer to the known p_ratio
#This is necessary because due to the nature of grid optimization, you must have discrete values that cover the range of parameter bounds
#Scipy minimize allows deviation from those discrete values and provides optimized parameter values over an entire range of bounds, not just the discrete values

#Almost identical to the earlier "prob_surv" function
#Instead of returning "None", it returns infinity for overflow errors
    #Because minimize optimizes the parameters based on the lowest value, this ensures that any parameter set returning infinity will not be considered optimal
def prob_surv_min(b, c, d, f, time):
    n_max = 100
    summation = 0
    for n in range(0, n_max):
        nfac = math.factorial(n)
        beta = (((-b) ** n) * (time ** ((c * n) + 1))) / ((c * n * nfac) + nfac)
        summation += beta
    try:
        survival_probability = math.exp(-d * time - f * summation)
    except OverflowError:
        survival_probability = float("inf") #Returning infinity

    return survival_probability

#Again, almost identical to the earlier "pratio_function" function
#Instead of returning "None", it returns infinity for overflow errors
#Additionally, instead of returning a p-ratio, it returns the absolute value of the difference between the p-ratio for a set of parameters and the known p-ratio
    #This value is ultimately the value that gets "minimized" and thus the parameters get optimized
def pratio_function_min(params):

    b_alt, c_alt, d_alt, f_alt, b_dos, c_dos, d_dos, d_non, f_non, alpha_alt, alpha_dos = params

    if alpha_alt + alpha_dos <= 1:
    
        f_dos = d_dos * -1
        alpha_non = 1 - alpha_alt - alpha_dos

        #Returning infinity for all overflow errors
        try:
            s_alt_time1 = prob_surv_min(b_alt,c_alt,d_alt,f_alt,time1)
        except OverflowError:
            s_alt_time1 = float("inf")

        try:
            s_alt_time2 = prob_surv_min(b_alt,c_alt,d_alt,f_alt,time2)
        except OverflowError:
            s_alt_time2 = float("inf")
            
        try:
            s_dos_time1 = prob_surv_min(b_dos,c_dos,d_dos,f_dos,time1)
        except OverflowError:
            s_dos_time1 = float("inf")

        try:
            s_dos_time2 = prob_surv_min(b_dos,c_dos,d_dos,f_dos,time2)
        except OverflowError:
            s_dos_time2 = float("inf")

        try:
            s_non_time1 = prob_surv_min(0,1,d_non,f_non,time1)
        except OverflowError:
            s_non_time1 = float("inf")

        try:
            s_non_time2 = prob_surv_min(0,1,d_non,f_non,time2)
        except OverflowError:
            s_non_time2 = float("inf")

        try: 
            temp = ((2*alpha_alt*s_alt_time1*s_alt_time2) + (2*alpha_dos*s_dos_time1*s_dos_time2) + (2*alpha_non*s_non_time1*s_non_time2)) / (((1-s_alt_time1)*alpha_alt*s_alt_time2) + ((1-s_dos_time1)*alpha_dos*s_dos_time2) + ((1-s_non_time1)*alpha_non*s_non_time2))
        except OverflowError:
            temp = float("inf")

        try:
            temp1 = (((1-s_alt_time1)*alpha_alt) + ((1-s_dos_time1)*alpha_dos) + ((1-s_non_time1)*alpha_non)) / ((2*alpha_alt*s_alt_time1) + (2*alpha_dos*s_dos_time1) + (2*alpha_non*s_non_time1))
        except OverflowError:
            temp1 = float("inf")

        try:
            pratio = temp * temp1
        except OverflowError:
            pratio = float("inf")

        #Returning the absolute value of difference between the p-ratio for a set of parameters and the known p-ratio
        return abs(pratio - known_p_ratio)

    else:
        return float("inf")


#Defining the bounds of the paramaters
#Bounds are RANGES of values instead of discrete values like earlier
#Additionally, instead of putting a "cap" on certain parameters (which is necessary for grid optimization), some parameters have no bounds based on Konrad et. Al
#Parameters for neo-functionalization / sub-functionalization
b_alt_range_min = (0.0001,None)
c_alt_range_min = (0.0001,None)
d_alt_range_min = (0.0001,None)
f_alt_range_min = (0.0001,None)

#Parameters for Dosage
b_dos_range_min = (None,-0.0001)
c_dos_range_min = (0.0001,0.9999)
d_dos_range_min =  (-0.1,-0.00001)

#Parameters for Non-functionalization       
d_non_range_min = (10,None) 
f_non_range_min = (None,None)

#Paramaters for the proportion of the initial genome for each category above
alpha_alt_range_min = (0.01,0.99)
alpha_dos_range_min = (0.01,0.99)

#Putting these bounds in a list
bounds = [b_alt_range_min, c_alt_range_min, d_alt_range_min, f_alt_range_min, b_dos_range_min, c_dos_range_min, d_dos_range_min, d_non_range_min, f_non_range_min, alpha_alt_range_min, alpha_dos_range_min]


#Initial guess for the optimization
#Using the values from the grid optimization as the initial guesses
initial_guess = [best_b_alt, best_c_alt, best_d_alt, best_f_alt, best_b_dos, best_c_dos, best_d_dos, best_d_non, best_f_non, best_alpha_alt, best_alpha_dos]

#Optimizing using SciPy minimize
#Nelder-Mead is a direct search optimization method that works well with a large number of parameters and bounds for those parameters
min_result = minimize(pratio_function_min, initial_guess, bounds = bounds, method = "Nelder-Mead")

#Creating a final list of optimized parameter values and adding in the parameter values that were not "bounded" (i.e., did not need to be 'optimized')
optimized_parameter_list_v1 = list(min_result.x)
calc_p_ratio_v1 = pratio_function(*optimized_parameter_list_v1, time1, time2)
optimized_parameter_list_v1.insert(7, optimized_parameter_list_v1[6] * -1) #dosage balance "f" is equal to negative dosage balance "d"
optimized_parameter_list_v1.insert(8, 0.0) #non-functionalization "b" is always equal to 0.0
optimized_parameter_list_v1.insert(9, 1.0) #non-functionalization "c" is always equal to 1.0
optimized_parameter_list_v1.append(1 - optimized_parameter_list_v1[12] - optimized_parameter_list_v1[13]) #non-functionalization genome percentage is 1 - altfunctionalization percentage - dosage balance percentage)


with open("optimization_output.txt", "a") as file:

    #Creating empty space
    file.write("\n" + "\n")
    file.write("Optimized parameters after conducting off-grid optimization using SciPy minimize Nelder-Mead method:" + "\n")
    
    #Writing each parameter and its optimized value
    file.write("b_alt = " + str(optimized_parameter_list_v1[0]) + "\n")
    file.write("c_alt = " + str(optimized_parameter_list_v1[1]) + "\n")
    file.write("d_alt = " + str(optimized_parameter_list_v1[2]) + "\n")
    file.write("f_alt = " + str(optimized_parameter_list_v1[3]) + "\n")
    file.write("b_dos = " + str(optimized_parameter_list_v1[4]) + "\n")
    file.write("c_dos = " + str(optimized_parameter_list_v1[5]) + "\n")
    file.write("d_dos = " + str(optimized_parameter_list_v1[6]) + "\n")
    file.write("f_dos = " + str(optimized_parameter_list_v1[7]) + "\n")
    file.write("b_non = " + str(optimized_parameter_list_v1[8]) + "\n")
    file.write("c_non = " + str(optimized_parameter_list_v1[9]) + "\n")
    file.write("d_non = " + str(optimized_parameter_list_v1[10]) + "\n")
    file.write("f_non = " + str(optimized_parameter_list_v1[11]) + "\n")
    file.write("alpha_alt = " + str(optimized_parameter_list_v1[12]) + "\n")
    file.write("alpha_dos = " + str(optimized_parameter_list_v1[13]) + "\n")
    file.write("alpha_non = " + str(optimized_parameter_list_v1[14]) + "\n")

    #Writing the resulting p-ratio achieved from these optimized values
    file.write("\n" + "Resulting p-ratio: " + str(calc_p_ratio_v1) + "\n")

    #Writing the percent error between the "optimized" p-value and the known p-value
    file.write("Percent error: " + str(((abs(calc_p_ratio_v1 - known_p_ratio))/known_p_ratio)*100))
