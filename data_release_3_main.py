from data_release_1 import *
import csv 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from scipy import stats
import statistics
import pandas as pd 
#importing  the necessary libraries for the analysis of the virus data 

#goal is to make a scatter plot showing day vs active infections 

with open(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#3.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    reader = csv.reader(f)
    headers = next(reader) # Get the first row to get the headers of the csv
    for h in headers:
        print(h) #show the headers of the csv

#Conner change this path to data release 3

#with open("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    #reader = csv.reader(f)
    #headers = next(reader) # Get the first row to get the headers of the csv
    #for h in headers:
       #print(h) #show the headers of the csv

Virus_count.instantiate_from_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#3.csv")
#Virus_count.instantiate_from_csv(r"/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv")#instantiating virus_count objects from the csv file using the class method instantiate_from_csv
#Conner change this path to data release 3


day = [] #creating a list to store the day of each virus_count object
active_reported_daily_cases = [] #creating a list to store the active reported daily cases of each virus_count object

for virus_count in Virus_count.all_virus_count: #iterating through all the virus_count objects to get the day and active reported daily cases of each virus_count day and store them in the respective lists
    day.append(float(getattr(virus_count, "day")))
    active_reported_daily_cases.append(float(getattr(virus_count, "active_reported_daily_cases")))  

# Convert lists to NumPy arrays
X = np.array(day)
y = np.array(active_reported_daily_cases)

# Define the exponential function to fit the data
def exponential_func(x, I0, r):
    return I0 * np.exp(r * x)

# Fit the exponential function to the data
popt, pcov = curve_fit(exponential_func, X, y, p0=(1, 0.1))  # Initial guess for parameters a and b
print("I0,:", popt[0], "r:", popt[1])

# Calulate R0 v1
R0_v1 = 1 + popt[1] * 9 # Assuming an infectious period of 9 days
print("R0_v1:", R0_v1)

#Calculate R0 v2
g = np.exp(popt[1])
R0_v2 = np.power(g, 9)
print("R0_v2:", R0_v2)

#Calculate R0 
R0 = (R0_v1 + R0_v2) / 2
print("R0:", R0)


# Plot
plt.scatter(X, y, color='blue', label='All days')
plt.plot(X, exponential_func(X, *popt), color='red', label='Exponential fit')
plt.ylabel('Active Reported Daily Cases')
plt.xlabel('Day')
plt.title('Scatter Plot of Day vs Active Reported Daily Cases')
plt.legend()
plt.show()

S0 = 17899
E0 = 0
I0 = 1
R0 = 0
day = 70
N = 17900
sigma = 0.205
beta = 0.277
gamma = 1/9
def seir(day,beta,sigma,gamma,S0,E0,I0,R0,N,h):
    S_list = []
    E_list = []
    I_list = []
    R_list = []
    
    S_list.append(S0)
    E_list.append(E0)
    I_list.append(I0)
    R_list.append(R0)

    for virus_day in range(day):
        h = 1
        # Euler derivatives
        dS = -beta * S_list[virus_day] * I_list[virus_day] / N
        dE = beta * S_list[virus_day] * I_list[virus_day] / N - sigma * E_list[virus_day]
        dI = sigma * E_list[virus_day] - gamma * I_list[virus_day]
        dR = gamma * I_list[virus_day]

        # Euler updates
        S_list.append(S_list[virus_day] + h * dS)
        E_list.append(E_list[virus_day] + h * dE)
        I_list.append(I_list[virus_day] + h * dI)
        R_list.append(R_list[virus_day] + h * dR)

    return S_list, E_list, I_list, R_list

#Conner change this path to data release 3 
#data = pd.read_csv("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv")

data = pd.read_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#3.csv")
data = data["active reported daily cases"].tolist()

def grid_search(day, N, S0, E0, I0, R0, data):
    # parameter grids
    beta_values = np.linspace(0.2, 0.7)
    sigma_values = np.linspace(0.1, 0.3)
    gamma_values = np.linspace(0.05, 0.25)

    # step size for Euler
    h = 1

    # storage for search results
    SSE_list = []
    beta_list = []
    sigma_list = []
    gamma_list = []


    best_SSE = float("inf")
    best_beta = None
    best_sigma = None
    best_gamma = None

    # iterate through parameter combinations
    for b in beta_values:
        for s in sigma_values:
            for g in gamma_values:
                # run Euler method for these parameters
                S_list, E_list, I_list, R_list = seir(day, b, s, g, S0, E0, I0, R0, N, h)

                model_I = I_list[: len(data)]

                # compute sum of squared errors
                errors = [(model_I[i] - data[i]) ** 2 for i in range(len(model_I))]
                sse = sum(errors)

                # record results
                SSE_list.append(sse)
                beta_list.append(b)
                sigma_list.append(s)
                gamma_list.append(g)

                # update best if this is lowest SSE so far
                if sse < best_SSE:
                    best_SSE = sse
                    best_beta = b
                    best_sigma = s
                    best_gamma = g
    
    return best_beta, best_sigma, best_gamma, best_SSE, SSE_list, beta_list, sigma_list, gamma_list

# perform grid search using defined parameters
best_beta, best_sigma, best_gamma, best_SSE, SSE_list, beta_list, sigma_list, gamma_list = \
    grid_search(day, N, S0, E0, I0, R0, data)

# print out best parameters and associated SSE
print(f"Best parameters from grid search:\n  beta = {best_beta}\n  sigma = {best_sigma}\n  gamma = {best_gamma}\n  SSE = {best_SSE}")

# generate model output with optimal parameters
S_list_opt, E_list_opt, I_list_opt, R_list_opt = seir(day, best_beta, best_sigma, best_gamma, S0, E0, I0, R0, N, h=1)

# plot the optimal SEIR curves
plt.plot(range(day+1), S_list_opt, label="Susceptible")
plt.plot(range(day+1), E_list_opt, label="Exposed")
plt.plot(range(day+1), I_list_opt, label="Infectious")
plt.plot(range(day+1), R_list_opt, label="Recovered")

plt.xlabel("Day")
plt.ylabel("Population")
plt.title("SEIR Model (best fit)")
plt.legend()

plt.show()

#peak analysis after extending days
extended_days = 200  # run model well past the observed 70 days
S_ext, E_ext, I_ext, R_ext = seir(extended_days, best_beta, best_sigma, best_gamma, S0, E0, I0, R0, N, h=1)

peak_I = max(I_ext)
peak_day = I_ext.index(peak_I)
print(f"Extended simulation peak infectious count: {peak_I:.2f} on day {peak_day}")