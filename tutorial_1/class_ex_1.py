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

def sell_item(item_name, quantity):
    if item_name not in inventory.keys():
        print("Item not found.")
    elif inventory[item_name]["quantity"] >= quantity:
        inventory[item_name]["quantity"] -=quantity
        print(f"Total Price: {inventory[item_name]["price"] * quantity}")
    else:
        print("Not enough stock.")
        
sell_item("banana", 2)