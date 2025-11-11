# -*- coding: utf-8 -*-
import re
import pandas as pd
# import os


# ex 1
def convert_case(text):
    upper = text.upper()
    lower = text.lower()
    return [upper, lower]

path = "/Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_2/class_ex_2.py"

def extract_file_info(path):
    splits = re.split(r"[./]", path)
    file_name = splits[len(splits)-2]
    extension = splits[len(splits)-1]
    return [file_name, extension]

convert_case("AsdfFBHdfdf")
extract_file_info(path)

# ex 2
# os.chdir("C:Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_2")

df = pd.read_csv("students.txt", delimiter=',')
print(df.head())
print(df.iloc[:5,])
df[["Name", "Grade"]]
print(df["Grade"].sum() / len(df.iloc[:]))

subject_math = df["Subject"] == "Math"
print(df.loc[subject_math])

city_new_york = df["City"] == "New York"
grade_above_85 = df["Grade"] > 85
print(df.loc[city_new_york & grade_above_85])

print(df.loc[1, "Grade"])
print(df.iloc[1, 4])

print(df.loc[df["Name"] == "Charlie"]["Grade"].values[0])

print(df.loc[df["Name"] == "Diana"]["Grade"])
df.loc[df["Name"] == "Diana", "Grade"] = 90
print(df.loc[df["Name"] == "Diana"]["Grade"])

best_on_top = df.sort_values("Grade", ascending=False)
best_on_top.iloc[0:3]

by_city_and_age = df.sort_values(by=["City", "Age"], ascending=False)
print(by_city_and_age)

by_city_and_age["Distinction"] = by_city_and_age["Grade"] >= 90
print(by_city_and_age)

by_city_and_age.to_csv("students_updated.csv", index=False)
