from data_release_1 import *
import csv 
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
from scipy import stats
import statistics
import pandas as pd 
#importing data_release_1 class and the necessary libraries for the analysis of the virus data 

#goal is to make a scatter plot showing day vs active infections 

#with open(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#1.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    #reader = csv.reader(f)
    #headers = next(reader) # Get the first row to get the headers of the csv
    #for h in headers:
     #   print(h) #show the headers of the csv

with open("/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#1.csv", newline="") as f: #opening the csv file to get the headers of the csv file
    reader = csv.reader(f)
    headers = next(reader) # Get the first row to get the headers of the csv
    for h in headers:
        print(h) #show the headers of the csv) as f: #opening the csv file to get the headers of the csv file

#Virus_count.instantiate_from_csv(r"C:\Users\dance\OneDrive - University of Virginia\Computational BME\Module-2-Epidemics-SIR-Modeling-Broderick_Looney\Data\mystery_virus_daily_active_counts_RELEASE#1.csv")
Virus_count.instantiate_from_csv(r"/Users/connerlooney/Documents/GitHub/Module-2-Epidemics-SIR-Modeling-Broderick_Looney/Data/mystery_virus_daily_active_counts_RELEASE#1.csv")#instantiating virus_count objects from the csv file using the class method instantiate_from_csv

day = [] #creating a list to store the day of each virus_count object
active_reported_daily_cases = [] #creating a list to store the active reported daily cases of each virus_count object

for virus_count in Virus_count.all_virus_count: #iterating through all the virus_count objects to get the day and active reported daily cases of each patient and store them in the respective lists
    day.append(float(getattr(virus_count, "day")))
    active_reported_daily_cases.append(float(getattr(virus_count, "active_reported_daily_cases")))  

X = np.array(day)  # Independent variable
y = np.array(active_reported_daily_cases)   # Dependent variable

#model = LinearRegression()
#model.fit(np.array(X).reshape(-1, 1), y)

#slope = model.coef_[0] # Get the slope (coefficient) of the linear regression model
#intercept = model.intercept_ # Get the intercept of the linear regression model 
#r2 = model.score(np.array(X).reshape(-1, 1), y) # Get the R-squared value of the linear regression model

#equation = f"y = {slope:.2f}x + {intercept:.2f}\nR^2 = {r2:.2f}" # Create the equation of the line in the form of y = mx + b
#plt.text((X.min() + X.max()) / 2, y.max(), equation, color='red', fontsize=12, ha='center', verticalalignment='top') # Add the equation of the line to the scatter plot

plt.scatter(X, y, color='blue') # Create a scatter plot of day vs active infections
#plt.plot(X, model.predict(np.array(X).reshape(-1, 1)), color='red') # Add the linear regression line to the scatter plot
plt.xlabel('Day') # Create a label for the x-axis
plt.ylabel('Active Reported Daily Cases') #create label for y axis
plt.title('Scatter Plot of Day vs Active Reported Daily Cases') #create title
plt.show()
