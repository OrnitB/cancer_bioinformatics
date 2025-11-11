# -*- coding: utf-8 -*-
import re

def convert_case(text):
    upper = text.upper()
    lower = text.lower()
    return [upper, lower]

def extract_file_info(path):
    splits = re.split(r"[./]", path)
    file_name = splits[len(splits)-2]
    extension = splits[len(splits)-1]
    return [file_name, extension]

convert_case("AsdfFBHdfdf")
extract_file_info("/Users/ornitb/Documents/Technion/research_and_studies/courses/ביואינפורמטיקה של סרטן/cancer_bioinformatics/tutorial_2/class_ex_2.py")