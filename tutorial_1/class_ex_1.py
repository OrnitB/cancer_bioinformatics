# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 19:52:59 2024

@author: Ornit Bhonkar
"""
from seaborn import load_dataset

iris_df = load_dataset('iris')
print(iris_df.head())

total_sepal_length = 0
virginica_count = 0

for row in iris_df.itertuples():
    if row.species == "virginica":
        virginica_count += 1
        total_sepal_length += row.sepal_length

if virginica_count > 0:
    print(f"Average sepal length for virginica species: {total_sepal_length / virginica_count:.3f}")
else:
    print("Average sepal length for virginica species: 0")



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
        
def add_item(item_name, price, quantity):
    if item_name in inventory:
        inventory[item_name]["quantity"] += quantity
        inventory[item_name]["price"] = price
    else:
        inventory[item_name] = {"quantity": quantity, "price": price}
        
        
sell_item("banana", 2)
print(inventory)
add_item("kiwi", 1.2, 4)
add_item("apple", 0.1, 2)
print(inventory)