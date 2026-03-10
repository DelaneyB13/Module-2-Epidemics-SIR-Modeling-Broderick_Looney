#AI usage statement: ChatGPT was used to help with the creation of Euler's method
from data_release_2 import *
import csv 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from scipy import stats
import statistics
import pandas as pd 
#importing data_release_1 class and the necessary libraries for the analysis of the virus data 

#goal is to make a scatter plot showing day vs active infections 

<<<<<<< HEAD:data_release_2_main.py
#with open(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#1.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    #reader = csv.reader(f)
    #headers = next(reader) # Get the first row to get the headers of the csv
    #for h in headers:
     #   print(h) #show the headers of the csv

with open("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv", newline="") as f: #opening the csv file to get the headers of the csv file
=======
with open(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#1.csv", newline="") as f: #opening the csv file to get the headers of the csv file
>>>>>>> cc5000e594f5026314694130f1a1f586c38e188:data_release_1_main.py
    reader = csv.reader(f)
    headers = next(reader) # Get the first row to get the headers of the csv
    for h in headers:
       print(h) #show the headers of the csv

<<<<<<< HEAD:data_release_2_main.py
#Virus_count.instantiate_from_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#2.csv")
Virus_count.instantiate_from_csv(r"/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv")#instantiating virus_count objects from the csv file using the class method instantiate_from_csv
=======
#with open("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#1.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    #reader = csv.reader(f)
    #headers = next(reader) # Get the first row to get the headers of the csv
    #for h in headers:
     #   print(h) #show the headers of the csv) as f: #opening the csv file to get the headers of the csv file

Virus_count.instantiate_from_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#1.csv")
#Virus_count.instantiate_from_csv(r"/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#1.csv")#instantiating virus_count objects from the csv file using the class method instantiate_from_csv
>>>>>>> 0cc5000e594f5026314694130f1a1f586c38e188:data_release_1_main.py

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

#Euler's method to find SIER model parameters
beta = 0.3     # transmission rate
sigma = 0.2    # incubation rate
gamma = 0.1    # recovery rate
h = 1          # time step

N = 1000000    # total population

S = N - 1
E = 0
I = 1
R = 0
seir_results = []

for virus_day in Virus_count.all_virus_count:


    day = virus_day.get_day()

    # Euler derivatives
    dS = -beta * S * I / N
    dE = beta * S * I / N - sigma * E
    dI = sigma * E - gamma * I
    dR = gamma * I

    # Euler updates
    S = S + h * dS
    E = E + h * dE
    I = I + h * dI
    R = R + h * dR

    seir_results.append((day, S, E, I, R))
    import matplotlib.pyplot as plt

S_list = []
E_list = []
I_list = []
R_list = []
days = []

for virus_day in Virus_count.all_virus_count:

    day = virus_day.get_day()

    dS = -beta * S * I / N
    dE = beta * S * I / N - sigma * E
    dI = sigma * E - gamma * I
    dR = gamma * I

    S = S + h*dS
    E = E + h*dE
    I = I + h*dI
    R = R + h*dR

    S_list.append(S)
    E_list.append(E)
    I_list.append(I)
    R_list.append(R)
    days.append(day)

plt.plot(days, S_list, label="Susceptible")
plt.plot(days, E_list, label="Exposed")
plt.plot(days, I_list, label="Infectious")
plt.plot(days, R_list, label="Recovered")

plt.xlabel("Day")
plt.ylabel("Population")
plt.title("SEIR Model")
plt.legend()

plt.show()