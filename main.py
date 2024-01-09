import subprocess
import time

file_path = './input_file_path.txt' 

with open(file_path, 'r', encoding='utf-8') as file:
    lines = file.readlines()

top_file_path = "./topol.top"
pdb_file_path = lines[6]
mdp_file_path = lines[9]
itp_files  = lines[12:]

itp_files_name = []

for itp_file in itp_files:
    itp_files_name.append(itp_file.split('/')[-1].strip())


input_molecule = pdb_file_path.replace('.pdb', '')
print(f"input_molecule : {input_molecule}") 
command = f'gmx pdb2gmx -f {input_molecule}.pdb -o {input_molecule}_processed.gro -water spce -ff oplsaa -p topol.top'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("pdb2gmx 실행 완료")
except subprocess.CalledProcessError as e:
    print(f"pdb2gmx 실행 중 오류 발생: {e}")

with open(top_file_path, 'r', encoding='utf-8') as top_file:
    top_content = top_file.read()

itp_files_section = '\n'.join(itp_files_name)

index_of_forcefield = top_content.find('; Include forcefield parameters')
if index_of_forcefield != -1:
    updated_top_content = (top_content[:index_of_forcefield + len('; Include forcefield parameters')] +
                           '\n' + '#include \"' + itp_files_section + '\"' + 
                           top_content[index_of_forcefield + len('; Include forcefield parameters'):])
else:
    updated_top_content = top_content + '\n' + itp_files_section

with open(top_file_path, 'w', encoding='utf-8') as top_file:
    top_file.write(updated_top_content)

command = f'gmx editconf -f {input_molecule}_processed.gro -o {input_molecule}_newbox.gro -c -d 1.0 -bt cubic'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("Box 생성 완료")
except subprocess.CalledProcessError as e:
    print(f"Box 생성 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx solvate -cp {input_molecule}_newbox.gro -cs spc216.gro -o {input_molecule}_solv.gro -p topol.top'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("물 분자 추가 완료")
except subprocess.CalledProcessError as e:
    print(f"물 분자 추가 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx grompp -f {mdp_file_path} -c {input_molecule}_solv.gro -p topol.top -o em.tpr'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("에너지 최적화 준비 완료")
except subprocess.CalledProcessError as e:
    print(f"에너지 최적화 준비 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx mdrun -v -deffnm em'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("에너지 최적화 완료")
except subprocess.CalledProcessError as e:
    print(f"에너지 최적화 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx grompp -f {mdp_file_path} -c em.gro -p topol.top -o nvt.tpr'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("NVT 준비 완료")
except subprocess.CalledProcessError as e:
    print(f"NVT 준비 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx mdrun -v -deffnm nvt'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("NVT 완료")
except subprocess.CalledProcessError as e:
    print(f"NVT 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx grompp -f {mdp_file_path} -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("NPT 준비 완료")
except subprocess.CalledProcessError as e:
    print(f"NPT 준비 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx mdrun -v -deffnm npt'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("NPT 완료")
except subprocess.CalledProcessError as e:
    print(f"NPT 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx grompp -f {mdp_file_path} -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("MD 준비 완료")
except subprocess.CalledProcessError as e:
    print(f"MD 준비 중 오류 발생: {e}")
time.sleep(1)

command = f'gmx mdrun -v -deffnm md_0_1'
print(f"Command : {command}")
try:
    subprocess.run(command, shell=True, check=True)
    print("MD 완료")
except subprocess.CalledProcessError as e:
    print(f"MD 중 오류 발생: {e}")
time.sleep(1)

user_input = input("계산할 것을 선택하세요. (1: 에너지, 2: RMSD, 3: RMSF, 4: Gyrate) : ")
input_list = user_input.split()

if  '1' in input_list:
    command = f'gmx energy -f md_0_1.edr -o potential.xvg'
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print("에너지 계산 완료")
    except subprocess.CalledProcessError as e:
        print(f"에너지 계산 중 오류 발생: {e}")
time.sleep(10)

if '2' in input_list:
    command = f'gmx rms -s md_0_1.tpr -f md_0_1.xtc -o rmsd.xvg'
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print("RMSD 계산 완료")
    except subprocess.CalledProcessError as e:
        print(f"RMSD 계산 중 오류 발생: {e}")
time.sleep(10)

if '3' in input_list:
    command = f'gmx rmsf -s md_0_1.tpr -f md_0_1.xtc -o rmsf.xvg'
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print("RMSF 계산 완료")
    except subprocess.CalledProcessError as e:
        print(f"RMSF 계산 중 오류 발생: {e}")
time.sleep(10)

if '4' in input_list:
    command = f'gmx gyrate -s md_0_1.tpr -f md_0_1.xtc -o gyrate.xvg'
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print("Gyrate 계산 완료")
    except subprocess.CalledProcessError as e:
        print(f"Gyrate 계산 중 오류 발생: {e}")
time.sleep(10)


subprocess.run('rm topol.top*', shell=True, check=True)
subprocess.run('rm mdout.mdp*', shell=True, check=True)
subprocess.run('rm posre.itp*', shell=True, check=True)
subprocess.run('rm conf.gro*', shell=True, check=True)
subprocess.run('rm out.gro*', shell=True, check=True)

print("Job Success!")