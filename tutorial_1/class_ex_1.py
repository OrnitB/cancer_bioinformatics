# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 19:52:59 2024

@author: Ornit Bhonkar
"""

inventory = {"apple": {"price": 0.5, "quantity": 10}, "banana": {"price": 0.3, "quantity": 15},\
             "chocolate": {"price": 1.5, "quantity": 5}}
    
print(inventory["apple"]["price"])
print(inventory["apple"]["quantity"])

def display_inventory():
    for key, val in inventory.items():
        print(f"{key} costs {val["price"]}, {val["quantity"]} left")
display_inventory()