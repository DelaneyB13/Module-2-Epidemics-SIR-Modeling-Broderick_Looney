from data_release_1 import *
import csv 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from scipy import stats
import statistics
import pandas as pd 
#importing  the necessary libraries for the analysis of the virus data 


#with open(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#2.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    #reader = csv.reader(f)
    #headers = next(reader) # Get the first row to get the headers of the csv
    #for h in headers:
       #print(h) #show the headers of the csv



with open("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    reader = csv.reader(f)
    headers = next(reader) # Get the first row to get the headers of the csv
    for h in headers:
       print(h) #show the headers of the csv

#Virus_count.instantiate_from_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#2.csv")
Virus_count.instantiate_from_csv(r"/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv")#instantiating virus_count objects from the csv file using the class method instantiate_from_csv



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


data = pd.read_csv("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv")

#data = pd.read_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#2.csv")
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

#Calculating Error 
# Et = true - approximate
true = 3294
approximate = peak_I
Et = true - approximate
relative_error = abs(Et) / abs(true) * 100
print(f"Relative error of the peak infectious count: {relative_error:.2f}%")

true_day = 83 
Et_day = true_day - peak_day
relative_error_day = abs(Et_day) / abs(true_day) * 100
print(f"Relative error of the peak day: {relative_error_day:.2f}%")
#end data release 2 

data_2 = pd.read_csv("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#2.csv")
#data_2 = pd.read_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#2.csv")
data_2 = data_2["active reported daily cases"].tolist()
data_3 = pd.read_csv("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#3.csv")
#data_3 = pd.read_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#3.csv")
data_3 = data_3["active reported daily cases"].tolist()
# compare the full release #3 dataset against the SEIR model using best parameters
model_S2, model_E2, model_I2, model_R2 = seir(121, best_beta, best_sigma, best_gamma, S0, E0, I0, R0, N, h=1)

x = np.arange(1, 123)  # Days 1 to 121
plt.figure()
plt.plot(x, model_I2, label="Model Infectious")
plt.scatter(range(len(data_2)), data_2, color='red', label="Observed active cases")
plt.xlabel("Day")
plt.ylabel("Infectious / Reported Active Cases")
plt.title("Predicted Infectious Model vs Data Release #2")
plt.axvline(x=73, color='gray', linestyle='--', label="Peak Infected Population (Day 73)")
plt.legend()
plt.show()


x = np.arange(1, 123)  # Days 1 to 121
plt.figure()
plt.plot(x, model_I2, label="Model Infectious")
plt.scatter(range(len(data_3)), data_3, color='green', label="Observed active cases (Release #3)")
plt.xlabel("Day")
plt.ylabel("Infectious / Reported Active Cases")
plt.title("Predicted Infectious Model vs Data Release #3")
plt.axvline(x=73, color='gray', linestyle='--', label="Peak Infected Population (Day 73)")
plt.legend()
plt.show()

#For recommending interventions at VT: VT student population as S0 = 31035, I0 = 1, R0 = 0, E0(VT) = E0(UVA), use same best fit parameters to model the epidemic for first 70 days 
S0_VT = 31035
E0_VT = E0
I0_VT = 1
R0_VT = 0

N_VT = S0_VT + E0_VT + I0_VT + R0_VT

model_S_VT, model_E_VT, model_I_VT, model_R_VT = seir(day, best_beta, best_sigma, best_gamma, S0_VT, E0_VT, I0_VT, R0_VT, N_VT, h=1)

x = np.arange(1, 72)  # Days 1 to 71
plt.figure()
plt.plot(x, model_I_VT, label="Model Infectious (VT)")
plt.xlabel("Day")
plt.ylabel("Infectious Population")
plt.title("Predicted Infectious Model for VT Student Population")
plt.legend()
plt.show()



#Intervention 1: Mask Mandate - reduce beta (Transmission) by 40% vs no intervention, keeping sigma and gamma the same. Model the epidemic for 120 days to see the long-term effects of the mask mandate on the infectious population at VT.
days_mask = 120
model_S_VT_long, model_E_VT_long, model_I_VT_long, model_R_VT_long = seir(days_mask, best_beta, best_sigma, best_gamma, S0_VT, E0_VT, I0_VT, R0_VT, N_VT, h=1)

day0 = 70
S_70 = model_S_VT_long[day0]
E_70 = model_E_VT_long[day0]
I_70 = model_I_VT_long[day0]
R_70 = model_R_VT_long[day0]


beta_mask = best_beta * 0.6
days_post70 = days_mask - day0

baseline_segment = model_I_VT_long[day0:days_mask+1]

vaccinated = min(2000, S_70)
effective_vaccinated = 0.9 * vaccinated

S_70_vax = S_70 - effective_vaccinated
E_70_vax = E_70
I_70_vax = I_70
R_70_vax = R_70 + effective_vaccinated

closure_days = 14
beta_closure = best_beta * 0.2

S_close, E_close, I_close, R_close = seir(
    closure_days,
    beta_closure,
    best_sigma,
    best_gamma,
    S_70,
    E_70,
    I_70,
    R_70,
    N_VT,
    h=1
)

S_84 = S_close[-1]
E_84 = E_close[-1]
I_84 = I_close[-1]
R_84 = R_close[-1]

days_after = days_mask - (day0 + closure_days)

S_after, E_after, I_after, R_after = seir(
    days_after,
    best_beta,
    best_sigma,
    best_gamma,
    S_84,
    E_84,
    I_84,
    R_84,
    N_VT,
    h=1
)
I_school = I_close + I_after[1:]
model_S_mask, model_E_mask, model_I_mask, model_R_mask = seir(days_post70, beta_mask, best_sigma, best_gamma, S_70, E_70, I_70, R_70, N_VT, h=1)

model_S_vax, model_E_vax, model_I_vax, model_R_vax = seir(
    days_post70,
    best_beta,        
    best_sigma,
    best_gamma,
    S_70_vax,
    E_70_vax,
    I_70_vax,
    R_70_vax,
    N_VT,
    h=1
)




x = np.arange(day0, days_mask + 1)  # Days 70 to 120
plt.figure()
plt.plot(x, model_I_vax, label="Vaccine Campaign (Day 70)")
plt.plot(x, baseline_segment, label="Baseline (no intervention)")
plt.plot(x, model_I_mask, label="Model Infectious with Mask Mandate")
plt.plot(x, I_school, label="2-Week School Closure")
plt.xlabel("Day")
plt.ylabel("Infectious Population")
plt.title("Interventions starting at Day 70")
plt.legend()
plt.show()



