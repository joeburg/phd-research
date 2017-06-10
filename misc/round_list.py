#program to round a list of numbers 

decimal_place = input("Give the decimal place that you want to round to (\
for tenths give '10', for hundredths give '100', etc.): ")

def round_list(bond_list, decimal_place):
    for i in range(len(bond_list)):
        round_up = math.ceil(bond_list[i]*decimal_place)/decimal_place
        round_down = math.floor(bond_list[i]*decimal_place)/decimal_place
        if abs(bond_list[i]-round_up) < abs(bond_list[i]-round_down):
            bond_list[i] = round_up
        else:
            bond_list[i] = round_down

    return bond_list 
