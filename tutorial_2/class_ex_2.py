# -*- coding: utf-8 -*-
import re
import pandas as pd


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
