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
    out_put_Atoms = []
    out_put_Bonds = []
    out_put_Angles = []
    out_put_Dihedrals = []

    #Atoms 
    for i in range(len(current_unit['atoms'])-1):
        out_put_Atoms.append(current_unit['atoms'][i])

    for j in range(len(add_unit['atoms'])):
        add_unit['atoms'][j][0] = int(add_unit['atoms'][j][0]) + int(current_unit['atoms'][-1][0]) - 1
        out_put_Atoms.append(add_unit['atoms'][j][0])

    #Bonds
    for i in range(len(current_unit['bonds'])-1):
        out_put_Bonds.append(current_unit["bonds"][i])
    
    for j in range(len(add_unit['bonds'])):
        add_unit[j][1] = int(add_unit['bonds'][j][1])
        out_put_Bonds.append(add_unit[j])

    #Angles
    for i in range(len(current_unit['angles'])-1):
        out_put_Angles.append(current_unit[i])
    for j in range(len(add_unit['angles'])):
        add_unit[j][2] = int(add_unit[j][2]) + len(current_unit) - 1
        out_put_Angles.append(add_unit[j])

    #Dihedrals
    for i in range(len(current_unit['dihedrals'])-1):
        out_put_Dihedrals.append(current_unit[i])



    return out_put_Atoms, out_put_Bonds, out_put_Angles, out_put_Dihedrals

# Example usage
file_path = './URE.itp'
topology_data = read_topology_file(file_path)


number_of_atoms = topology_data.get('atoms')[-1][0]
n_times = 3



# Access the data using keys
print("Moleculetype:", topology_data.get('moleculetype'))
print("Atoms:", topology_data.get('atoms'))
print("Bonds:", topology_data.get('bonds'))
print("Angles:", topology_data.get('angles'))
print("Dihedrals:", topology_data.get('dihedrals'))

for i in range(n_times):
    out_put = add_unit(topology_data, topology_data)
    
