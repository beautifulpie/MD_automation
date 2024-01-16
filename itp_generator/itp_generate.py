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

def angle_parameter(i, j, k, atoms):
    # i, j, k : index of atoms
    # atoms : list of atoms
    theta = 0
    cth = 0
    return theta, cth

def dihedral_parameter(i, j, k, l, atoms):
    # i, j, k, l : index of atoms
    # atoms : list of atoms
    phi_0 = 0
    cp = 0
    mult = 0
    return phi_0, cp, mult

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
        inst[0][0] = int(inst[0][0]) 
        inst[0][2] = int(inst[0][2])
        inst[0][5] = int(inst[0][5])
        inst[0][6] = float(inst[0][6]) 
        inst[0][7] = float(inst[0][7])

        output_Atoms.append([int(inst[0][0])+ length_of_current - 1, inst[0][1], int(inst[0][2]), inst[0][3], inst[0][4], int(inst[0][5]) + length_of_current - 1, float(inst[0][6]), float(inst[0][7])])
        
    #Bonds
    locb = len(current_unit['bonds'])
    for i in range(locb):
        output_Bonds.append(current_unit["bonds"][i])
    
    for j in range(len(add_unit['bonds'])):
        add_unit['bonds'][j][0] = int(add_unit['bonds'][j][0])
        add_unit['bonds'][j][1] = int(add_unit['bonds'][j][1])
        output_Bonds.append([int(add_unit['bonds'][j][0]) + locb   , int(add_unit['bonds'][j][1]) + locb  ])

    #Angles
        
    for i in range(len(current_unit['angles'])-1):
        output_Angles.append(current_unit['angles'][i])

    output_Angles.append([ current_unit["angles"][-1][0], current_unit["angles"][-1][1] , len(current_unit["angles"]) + 1 , 
                        angle_parameter(current_unit["angles"][-1][0], current_unit["angles"][-1][1], current_unit["angles"][-1][2], [current_unit["atoms"][int(current_unit["angles"][-1][0])][1], current_unit["atoms"][int(current_unit["angles"][-1][1])][1], current_unit["atoms"][ int(current_unit["angles"][-1][2])][1]])[0], 
                        angle_parameter(current_unit["angles"][-1][0], current_unit["angles"][-1][1], current_unit["angles"][-1][2], [current_unit["atoms"][int(current_unit["angles"][-1][0])][1], current_unit["atoms"][int(current_unit["angles"][-1][1])][1], current_unit["atoms"][ int(current_unit["angles"][-1][2])][1]])[1]  
                        ])
    
    for j in range(len(add_unit['angles'])):
        output_Angles.append([int(add_unit['angles'][j][0]) + len(current_unit['angles']) - 1 , int(add_unit['angles'][j][1]) + len(current_unit['angles']), add_unit['angles'][j][2] , add_unit['angles'][j][3], add_unit['angles'][j][4]])

    #Dihedrals
    for i in range(len(current_unit['dihedrals'])):
        output_Dihedrals.append(current_unit['dihedrals'][i])

    return output_Atoms, output_Bonds, output_Angles, output_Dihedrals

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
print(f"Atoms : {out_put[0]}")
print(f"Bonds : {out_put[1]}")
print(f"Angles : {out_put[2]}")
print(f"Dihedrals : {out_put[3]}")