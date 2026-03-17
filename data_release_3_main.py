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
true = 3294 #found using the data from release #3, the true peak infectious count is 3294
approximate = peak_I #the peak infectious count from the extended SEIR model simulation
Et = true - approximate
relative_error = abs(Et) / abs(true) * 100 #formula given for relative error: |Et| / |true| * 100
print(f"Relative error of the peak infectious count: {relative_error:.2f}%")

true_day = 83 #found using the data from release #3, the true peak day is day 83
Et_day = true_day - peak_day
relative_error_day = abs(Et_day) / abs(true_day) * 100 #formula given for relative error: |Et| / |true| * 100
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

#Plotting the model predictions against the data from release #2 and release #3, showing the peak day and peak infectious population from the model on the plot for reference across the full 120 days of the epidemic, and comparing the model predictions to the observed data from both releases to evaluate how well the model captures the epidemic dynamics over time, including the timing and magnitude of the peak infectious population.
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
#the values were given from the powerpoint in class 

N_VT = S0_VT + E0_VT + I0_VT + R0_VT #changing N for model to reflect the VT student population instead of the UVA student population

model_S_VT, model_E_VT, model_I_VT, model_R_VT = seir(day, best_beta, best_sigma, best_gamma, S0_VT, E0_VT, I0_VT, R0_VT, N_VT, h=1) #running the SEIR model for the VT student population using the best fit parameters from the grid search

#using the model predictions for the VT student population to plot the predicted infectious population over the first 70 days of the epidemic, which can be used to inform intervention strategies at VT by showing the expected trajectory of the epidemic in the absence of interventions and identifying when the peak infectious population is expected to occur, which can help guide timing of interventions to mitigate the spread of the virus among the VT student population.
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

#making sure that the model are going off of the same initial conditions at day 70 for all interventions to ensure a fair comparison of the interventions starting at the same point in the epidemic, which allows us to evaluate the relative effectiveness of each intervention in reducing the infectious population at VT over time.
day0 = 70
S_70 = model_S_VT_long[day0]
E_70 = model_E_VT_long[day0]
I_70 = model_I_VT_long[day0]
R_70 = model_R_VT_long[day0]

#change to beta from the mask mandate intervention, which reduces transmission by 40%, so beta is reduced to 60% of the original best beta value from the grid search, and then model the epidemic for the remaining days after day 70 to see how the mask mandate affects the trajectory of the infectious population at VT compared to the baseline scenario with no intervention.
beta_mask = best_beta * 0.6
days_post70 = days_mask - day0


baseline_segment = model_I_VT_long[day0:days_mask+1] #the baseline segment of the infectious population trajectory from day 70 to day 120 without any interventions, which serves as a reference point for comparing the effects of the mask mandate and other interventions on the infectious population at VT over time.


#modeling the epidemic for the mask mandate intervention, where we use the adjusted beta value that reflects the 40% reduction in transmission due to the mask mandate, and the same sigma and gamma values from the best fit parameters, starting from the same initial conditions at day 70, to see how the mask mandate affects the trajectory of the infectious population at VT compared to the baseline scenario with no intervention and the other interventions.
model_S_mask, model_E_mask, model_I_mask, model_R_mask = seir(days_post70, beta_mask, best_sigma, best_gamma, S_70, E_70, I_70, R_70, N_VT, h=1)



#intervention 2: Vaccine Campaign - vaccinate 2000 students at day 70, with 90% vaccine efficacy, which effectively reduces the susceptible population by the number of effective vaccinations, and increases the recovered population by the same amount, while keeping the exposed and infectious populations the same at day 70. Then model the epidemic for the remaining days after day 70 to see how the vaccine campaign affects the trajectory of the infectious population at VT compared to the baseline scenario with no intervention and the mask mandate intervention.
vaccinated = min(2000, S_70) #the number of students vaccinated at day 70, which is the minimum of 2000 or the current susceptible population at day 70 to ensure that we do not vaccinate more students than are currently susceptible in the model
effective_vaccinated = 0.9 * vaccinated #the number of effective vaccinations, which is the number of vaccinated students multiplied by the vaccine efficacy of 90%, representing the reduction in the susceptible population due to the vaccine campaign, and the corresponding increase in the recovered population at day 70 in the SEIR model.


#adjusting the susceptible and recovered populations at day 70 to reflect the effects of the vaccine campaign, where the susceptible population is reduced by the number of effective vaccinations, and the recovered population is increased by the same amount, while keeping the exposed and infectious populations the same at day 70 to model the impact of the vaccine campaign on the epidemic trajectory at VT starting from day 70.
S_70_vax = S_70 - effective_vaccinated
E_70_vax = E_70
I_70_vax = I_70
R_70_vax = R_70 + effective_vaccinated

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


#intervention 3: School Closure - close schools for 2 weeks starting at day 70, which reduces the transmission rate (beta) by 80% during the closure period, and then model the epidemic for the remaining days after day 70 to see how the school closure affects the trajectory of the infectious population at VT compared to the baseline scenario with no intervention, the mask mandate intervention, and the vaccine campaign intervention.
closure_days = 14 #the duration of the school closure intervention, which is 14 days (2 weeks) starting from day 70, during which the transmission rate (beta) is reduced by 80% to model the impact of the school closure on the epidemic trajectory at VT during that period, and then returns to the original beta value after the closure period ends to model the subsequent trajectory of the infectious population at VT after the schools reopen.
beta_closure = best_beta * 0.2 #the adjusted beta value during the school closure period, which is 20% of the original best beta value from the grid search to reflect the 80% reduction in transmission due to the school closure, and then model the epidemic for the closure period and the subsequent days after day 70 to see how the school closure affects the trajectory of the infectious population at VT compared to the other interventions and the baseline scenario with no intervention.


#modeling the epidemic for the school closure intervention, where we first model the epidemic for the closure period with the reduced beta value to see how the school closure affects the infectious population at VT during that period, and then model the epidemic for the remaining days after day 70 with the original best beta value to see how the infectious population at VT evolves after the schools reopen, allowing us to evaluate the overall impact of the school closure intervention on the epidemic trajectory at VT compared to the other interventions and the baseline scenario with no intervention.
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

#after modeling the epidemic for the school closure intervention, we take the final values of the susceptible, exposed, infectious, and recovered populations at the end of the closure period (day 84) as the new initial conditions for modeling the epidemic for the remaining days after day 70 with the original best beta value to see how the infectious population at VT evolves after the schools reopen, allowing us to evaluate the overall impact of the school closure intervention on the epidemic trajectory at VT compared to the other interventions and the baseline scenario with no intervention.
S_84 = S_close[-1]
E_84 = E_close[-1]
I_84 = I_close[-1]
R_84 = R_close[-1]

#modeling the epidemic for the remaining days after day 70 with the original best beta value after the school closure period ends, using the final values from the school closure model as the new initial conditions, to see how the infectious population at VT evolves after the schools reopen, allowing us to evaluate the overall impact of the school closure intervention on the epidemic trajectory at VT compared to the other interventions and the baseline scenario with no intervention.
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

I_school = I_close + I_after[1:] #the total infectious population trajectory for the school closure intervention, which combines the infectious population during the closure period (I_close) with the infectious population after the closure period ends (I_after) to show the overall impact of the school closure intervention on the infectious population at VT over time compared to the other interventions and the baseline scenario with no intervention.




#Finally, we can plot the infectious population trajectories for the baseline scenario with no intervention, the mask mandate intervention, the vaccine campaign intervention, and the school closure intervention starting from day 70 to day 120 to visually compare the effects of each intervention on the trajectory of the infectious population at VT over time, which can help inform decision-making about which interventions may be most effective in mitigating the spread of the virus among the VT student population.

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



