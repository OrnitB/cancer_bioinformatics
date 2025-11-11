# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 19:52:59 2024

@author: stavn
"""

##Conditions
#If-else
x=80
y=40
if x==y*2:
    print("True")
else:
    print("False")
print(x)
#If-else-if
x=5
y=4
z=3

if x==y:
    print("x equals y")
else:
    if x==z:
        print("x equals z")
    else:
        if y==z:
            print("y equals z")
        else:
            print("all elements are different")
#elif
if x==y:
    print("x equals y")
elif x==z:
    print("x equals z")
elif y==z:
    print("y equals z")
else:
    print("all elements are different")
    
#match-case vs elif
x=5
if x==1:
    print("x=1")
elif x==2:
    print("x=2")
elif x==3:
    print("x=3")
elif x==4:
    print("x=4")
elif x==5:
    print("x=5")
else:
    print("x is not between 1-5")
    
match x:
    case 1:
        print("x=1")
    case 2:
        print("x=2")
    case 3:
        print("x=3")
    case 4:
        print("x=4")
    case 5:
        print("x=5")
    case _:
        print("x is not between 1-5")
            
#Ternary operator vs if-else
answer = True if x>=1 else False

if x>=1:
    answer = True
else:
    answer = False
##Loops
#While loop
x = 500
while x != 0:
    print("x = " + str(x))
    x//=2
#What's the last printed X? 
#What's x value?

#break

i = 0
while i < 9:
  i += 1
  if i % 3==0:
    break
  print(i)
  
#continue

i = 0
while i < 9:
  i += 1
  if i == 3:
    break
  print(i)

#for loop- using range() function only with 'stop' argument
for i in range(6):
    print(i)
#using range function with start,stop,step arguments
for i in range(6,2,-1):
    print(i)
#iterating over a list with enumerate
my_list = ['apple', 'banana', 'cherry', 'date']
for i,item in enumerate(my_list):
    print("i = "+ str(i) + ", item= " + item)
    

  
my_list = ['apple', 'banana', 'cherry', 'date']
print(my_list[0])
print(my_list[1])


my_list = ['apple', 'banana', 'cherry', 'date','watermelon','peach']
print(my_list[:3])
print(my_list[2:5])

print(list(enumerate(my_list)))

##Functions
#A function without arguments
def eize_kef_li():
    print("Woohoo!")

eize_kef_li()
#A function with arguments
def eize_kef_to(name):
    print("My name is " + name + ", eize kef li! Woohoo!")

eize_kef_to("Stav")
#A function that returns a value
def func1(arg1,arg2): 
    return arg1 - arg2 / 2    
func1(5, 3)
#Call the function
var1 = 5
var2 = 10
result = func1(arg1 = var1, arg2 = var2)
result
result = func1(var1, var2)
result


###Data structures
##Lists

my_list= ['one',2,'three','four','five']
print(my_list)

#reassign an element
my_list[1]='two'
print(my_list)

#append- adds a single element at the end of a list
my_list.append('six')
print(my_list)

#insert- adds a single element at the specified position

my_list.insert(4,"4.5")
print(my_list)

#extend- concatenate lists
my_second_list = ['six','seven','eight']
my_list.extend(my_second_list)
print(my_list)

#count- returns the number of elements in list with specified value
print(my_list.count('six'))

#index- returns the index of the first element with specified value
print(my_list.index('six'))

#pop- removes the element at the specified position
#also returns removed value
my_list.pop(4)
print(my_list)

#remove- removes the first item with the specified value
my_list.remove("six")
print(my_list)

#reverse- return the list in reverse order
my_list.reverse()
print(my_list)

#sort- sort the list - DOESN'T ALWAYS WORK
my_list.sort(reverse=True)
print(my_list)

##tuples
my_tuple = ('one',2,'three','four','five')
print(my_tuple)
my_tuple[1] = 'two'
print(my_tuple)



##sets
my_set = {1,1,3,4,5}
print(my_set)
#Trying to get to a set item- fails
my_set[1]
my_set.update("5")
print(my_set)
#You cannot access items in a set through index, because it
#has no order.
#But, you can loop through the items
for item in my_set:
    print(item)
my_set_2 = {1,3,6,7}
#set operations
print(my_set.union(my_set_2)) # איחוד בין 2 הקבוצות, כל האיברים כולל החיתוך שלהם
print(my_set.intersection(my_set_2)) # רק החיתוך בין הקבוצות
print(my_set.symmetric_difference(my_set_2)) # הכל חוץ מהחיתוך
print(my_set.difference(my_set_2)) # חיסור, האיברים של הראשון ללא האיברים שגם יש בשני

##dicts
my_dict = {"Shas":"Haredim","Yehadut_Hatorah":"Haredim",\
           "Meretz":"Smolanim"}

my_dict.get("Shas")
#Iterate over dictionary items

for key,val in my_dict.items():
    print("key = " + key + ", val = " + val)
    

#Dictionary assignment
my_dict["Shas"] = "Mizrahim"
my_dict["Ha'avoda"] = "Smolanim"

for key,val in my_dict.items():
    print(key + " - " + val)

my_dict.keys()
my_dict.values()

##Functions



# Exercise 1 solution
inventory = {
    "apple": {"price": 0.5, "quantity": 10},
    "banana": {"price": 0.3, "quantity": 15},
    "chocolate": {"price": 1.5, "quantity": 5}
}

inventory["apple"]["price"]

def display_inventory():
    for item in inventory:
        elem = inventory.get(item)
        print(item + "-" + str(elem["price"]) + "$, Quantity: "\
              + str(elem["quantity"]) )

def sell_item(item_name,quantity):
    if item_name not in inventory.keys():
        print("Error- Item does not exist")
        return
    if quantity <= inventory[item_name]["quantity"] :
        inventory[item_name]["quantity"] -=quantity
    else:
        print("Not enough stack")
        
def add_item(item_name,price,quantity):
    if item_name in inventory.keys():
        inventory[item_name]["price"] = price
        inventory[item_name]["quantity"] += quantity
    else:
        inventory[item_name] = {"price":price, "quantity":quantity}

    


