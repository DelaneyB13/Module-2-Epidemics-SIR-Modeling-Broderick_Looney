import csv # Importing the csv to get data from the CSV file

class Virus_count : #Create a class called virus_count
    all_virus_count = [] #create an empty list to store all virus_count instances

    def __init__(self, day: int, date: str, active_reported_daily_cases : int): #create a constructor that lists the attributes of the virus_count class needed for our analysis
        
        self.day = day  
        self.date = date
        self.active_reported_daily_cases = active_reported_daily_cases

        Virus_count.all_virus_count.append(self) #append the virus_count instance to the class to update the class variable all_virus_count with the new virus_count instance created   

    def __repr__(self): #creating a representer that defines what is shown to us when we print a virus_count type object
        return f"( {self.day} | {self.date} | {self.active_reported_daily_cases})" #this contains the 3 categories needed to complete the anaylsis of the virus_count data
    
    
    def get_day(self): # creating a getter method to get the day of a virus_count object
        return self.day

    def get_active_reported_daily_cases(self): # creating a getter method to get the active reported daily cases of a virus_count object
        return self.active_reported_daily_cases

    
    @classmethod #creating a class method to instantiate virus_count objects from a csv file
    def instantiate_from_csv(cls, filename: str):
        with open(filename, encoding="utf8") as f:
            reader = csv.DictReader(f)
            rows_of_virus_days = list(reader)
        for row in rows_of_virus_days:
            #creating a try except block to fix any value errors that may occur when trying to create virus_count objects from the csv file, if there is a value error it will just skip that row and move on to the next one
            try:    
                Virus_count(
                day = int(row['day']),
                date = row['date'],
                active_reported_daily_cases = int(row['active reported daily cases'])
            )
            except ValueError:
                continue

    @classmethod #creating a class method to get a virus_count object based on the day of the virus_count (what we are interested in)
    def get_virus_count_by_day(cls, day):
        for virus_count in Virus_count.all_virus_count:
            if day == virus_count.day: 
                return virus_count
            

    @classmethod #creating a class method to get a virus_count object based on the active reported daily cases of the virus_count (what we are interested in)
    def get_virus_count_by_active_cases(cls, active_reported_daily_cases):
        for virus_count in Virus_count.all_virus_count:
            if active_reported_daily_cases == virus_count.active_reported_daily_cases:
                return virus_count
            
    @classmethod #creating a class method filter to filter virus_count objects based on the attributes of the virus_count class
    def filter(cls, day: int = None, active_reported_daily_cases: int = None):
        all_virus_count = cls.all_virus_count
        remove_list = []
        attr_list = (
                        day,    
                        active_reported_daily_cases,
                        )
        attr_name = (
                        "day",
                        "active_reported_daily_cases",
                        )
        for attr in range(len(attr_list)):
            if attr_list[attr] is not None:
                for virus_count in all_virus_count:
                    if getattr(virus_count,attr_name[attr]) != attr_list[attr]:
                        remove_list.append(virus_count)
                all_virus_count = [virus_count for virus_count in all_virus_count if virus_count not in remove_list]
                remove_list.clear()
        return all_virus_count
    
        