def is_numeric(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def read_topology_file(file_path):
    data = {}
    current_key = None

    with open(file_path, 'r', encoding='UTF-8' ) as file:
        for line in file:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Check for section headers
            if line.startswith('[') and line.endswith(']'):
                current_key = line[2:-2].lower()
                data[current_key] = []
            elif line.startswith(';'):
                continue
            elif line.startswith('#'):
                continue
            else:
                # Split the line into tokens
                tokens = line.split()
                semicolon_index = tokens.index(';') if ';' in tokens else None
                if semicolon_index is not None:
                    tokens = tokens[:semicolon_index]
                data[current_key].append(tokens)

    return data

def add_unit(current_unit, add_unit ):
    output_Atoms = []
    output_Bonds = []
    output_Angles = []
    output_Dihedrals = []
    length_of_current = int(current_unit['atoms'][-1][0])
    #Atoms 
    for i in range(length_of_current-1):
        output_Atoms.append(current_unit['atoms'][i])
    
    for j in range(int(add_unit['atoms'][-1][0])):
        inst = []
        inst.append(add_unit['atoms'][j])
        #inst[0][0] = str(int(inst[0][0]) + length_of_current - 1)
        output_Atoms.append([int(inst[0][0])+ length_of_current - 1, inst[0][1], int(inst[0][2]), inst[0][3], inst[0][4], int(inst[0][5]) + length_of_current - 1, float(inst[0][6]), float(inst[0][7])])
        
    #Bonds
    for i in range(len(current_unit['bonds'])-1):
        output_Bonds.append(current_unit["bonds"][i])
    
    for j in range(len(add_unit['bonds'])):
        add_unit['bonds'][j][1] = int(add_unit['bonds'][j][1]) + 1
        output_Bonds.append(add_unit['bonds'][j])

    #Angles
    for i in range(len(current_unit['angles'])-1):
        output_Angles.append(current_unit['angles'][i])
    for j in range(len(add_unit['angles'])):
        add_unit['angles'][j][2] = int(add_unit['angles'][j][2]) + len(current_unit) - 1
        output_Angles.append(add_unit['angles'][j])

    #Dihedrals
    #for i in range(len(current_unit['dihedrals'])-1):
       # out_put_Dihedrals.append(current_unit['diherals'][i])

    return output_Atoms, output_Bonds, output_Angles#, out_put_Dihedrals

# Example usage
file_path = './URE.itp'
topology_data = read_topology_file(file_path)


number_of_atoms = topology_data.get('atoms')[-1][0]
n_times = 1



# Access the data using keys
print("Moleculetype:", topology_data.get('moleculetype'))
print("Atoms:", topology_data.get('atoms'))
print("Bonds:", topology_data.get('bonds'))
print("Angles:", topology_data.get('angles'))
print("Dihedrals:", topology_data.get('dihedrals'))

for i in range(n_times):
    out_put = add_unit(topology_data, topology_data)


print()
print()
print()
print(out_put[0])
