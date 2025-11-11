
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 11:15:55 2024

@author: stavn
"""


import pandas as pd
import os

# os.chdir("C:Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_2")

# os.chdir(r"C:\Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_2")

os.getcwd()
##File handling

# Python code to create a file
file = open('Tuto1.txt','w')
file.write("I love Yosi's course! \n")
file.write("Stav is also OK..!")
file.close()

file = open('Tuto1.txt','a')
file.write("\nBut what about my grades?")
file.close()

file = open('Tuto1.txt','r')
for line in file.readlines():
    print(line)


##String manipulations
#Endswith, startswith
string = "Donald Trump"
string.startswith("Do")
string.endswith("mp")
#Isalpha, isalnum, isdigit
string = "asdfhabdksdf892341341234!@#$!@#$"
string.isalpha() 
"ACGTGTGCGAGCAGCA".isalpha()

"3ACTGAGTCGAGCAGT5".isalnum()

"5123451245ADFADF".isdigit()
#Upper, lower
string = "don't worry about a thing, every little thing is gonna be alright"
string.upper()
string = "DON'T STOP ME NOW!"
string.lower()

#strip
string = "/path/"
string.strip("/")
string = "____main____"
string.strip("_")

#split
string= "/path/to/file.txt"
string.split(sep="/")

#re.split- for multiple delimiters- in this case . and /
import re
re.split(r"[./]",string)

string = string.replace("path", "road")
print(string)

#Initialize data as dictionary
data = {"animal":["Dog","Cat","Donkey","Cow"],\
        "avg_weight_kg":[10,4,120,250],\
        "sound":["woof-woof","meow","eeh-aah","moo"],\
        "is_cute":[True,True,False,True]
        }
#Transform to data frame
df = pd.DataFrame(data)

#Access to rows via slicing
df[1:3]
#Access to a single row via slicing - FAILS
df[2]
#Access to rows via iloc() function
df.iloc[3]
df.iloc[1:3]

#Adding index names by specific column
df = df.set_index("animal")

#Access to rows via loc() function
df.loc["Dog"]
rows = ["Cat","Donkey"]
df.loc[rows]

#Conditional selection- Only the cute animals
df.loc[df["is_cute"] == True]

#Only the non-cute animals that weigh over 60 kg on average
df.loc[(df["avg_weight_kg"] >= 60) & (df["is_cute"] == False)]

##Access to columns

#Access to columns via column name
df["is_cute"]
cols = ["sound","avg_weight_kg"]
df[cols]

#Access to columns via attribute
df.sound
df.avg_weight_kg


#Access to specific elements
print(df.loc["Donkey","sound"])
rows = ["Donkey","Cow"]
cols = ["avg_weight_kg","is_cute"]
print(df.loc[rows,cols])
#Access to elements with slicing
features_list = ["sound","is_cute"]
print(df.loc["Cat",features_list])

#Combination- Use boolean indexing with column selection
#Filter only the animals whose sound starts with a letter bigger than "n"
#Filter only columns is_cute, sound, in that order
sounds_boolean = df["sound"]>"n"
cols = ["is_cute","sound"]
df.loc[sounds_boolean, cols]



#Exercise 1 String Manipulation solution:
#convert case
def convert_case(text):
    upper = text.upper()
    lower = text.lower()
    return [upper,lower]

#extract_file_info
path = "C:/Users/stavn/OneDrive/שולחן העבודה/Cancer_Bioinformatics/2024-2025/2_Table_manipulation/sort_by_grade_students.csv"

def extract_file_info(path):
    path_split = re.split(r"[./]",path)
    file_name = path_split[len(path_split)-2]
    extension = path_split[len(path_split)-1]
    return [file_name,extension]




#Exercise 2 Data Frames solution:
students = pd.read_csv("students.txt")
#Print the first five rows
print(students.loc[:4,])

#Print the average grade of all students
print(students[["Name","Grade"]])
print(students["Grade"].mean())

#Display students from New York who scored more than 85
print(students.loc[(students["City"] == "New York") & (students["Grade"] > 85)])

#Display all rows where subject is Math
print(students.loc[students["Subject"] == "Math",:])

#Print the grade of Charlie
print(students.loc[students["Name"] == "Charlie",:])

#Update the grade of Diana to 90 and print the updated row
students.loc[students["Name"] == "Diana","Grade"] = 90
print(students.loc[students["Name"] == "Diana",:])

#Sort the data frame in descending order and display top 3 rows
sort_by_grade_students = students.sort_values(by = "Grade",ascending = False)
sort_by_grade_students.iloc[0:3]
#Sort the DataFrame by City and then by Age within each city.
sort_by_city_students = students.sort_values(by = "City")
sort_by_age_and_city_students = students.sort_values(by=["City","Age"])

#Add a new column, Distinction, based on whether the grade is 90 or higher.
sort_by_grade_students["Distinction"] = sort_by_grade_students["Grade"] >= 90
#Save the modified data frame to a new file without indexing
sort_by_grade_students.to_csv("sort_by_grade_students.csv", index=True)
