import subprocess
import time
import md_automation as md

def main():
    file_path = './input_file_path.txt'
    lines = md.read_lines_from_file(file_path)

    top_file_path = "./topol.top"
    pdb_file_path = lines[6].strip()
    mdp_file_path = lines[9].strip('\n')
    minimization_mdp = lines[10].strip('\n')
    nvt_mdp = lines[11].strip('\n')
    npt_mdp = lines[12].strip('\n')
    itp_files = lines[15:]
    itp_files_name = md.extract_itp_files_names(itp_files)

    print("====================================")
    print('\033[1m'+"MD_automation"+'\0333')
    print("Made by Youngwoo Jung")
    print("Final update : 2024.01.09")
    print("====================================")

    print(f"Load pdb file : {pdb_file_path}")
    print(f"Load mdp file : {mdp_file_path}")
    print(f"Load minimization mdp file : {minimization_mdp}")
    print(f"Load nvt mdp file : {nvt_mdp}")
    print(f"Load npt mdp file : {npt_mdp}")

    print("====================================")

    input_molecule = pdb_file_path.replace('.pdb', '')
    print(f"input_molecule : {input_molecule}")

    # pdb2gmx 실행
    pdb2gmx_command = f'gmx pdb2gmx -f {input_molecule}.pdb -o {input_molecule}_processed.gro -water spce -p topol.top'
    subprocess.run(pdb2gmx_command, shell=True, check=True)
    print("pdb2gmx 완료")
    
    with open(top_file_path, 'r', encoding='utf-8') as top_file:
        top_content = top_file.read()

    itp_files_section = '\n'.join(itp_files_name)

    commands = [
        f'gmx editconf -f {input_molecule}_processed.gro -o {input_molecule}_newbox.gro -c -d 1.0 -bt cubic', # 박스 생성
        f'gmx solvate -cp {input_molecule}_newbox.gro -cs spc216.gro -o {input_molecule}_solv.gro -p topol.top', # 수분 분자 추가
    ]
    for command in commands:
        md.run_gmx_command(command, "작업 완료")
        md.print_rotating_bar()

#    index_of_forcefield = top_content.find('; Include forcefield parameters')
#    if index_of_forcefield != -1:
#        updated_top_content = (top_content[:index_of_forcefield + len('; Include forcefield parameters')] +
#                               '\n' + '#include \"' + itp_files_section + '\"' + 
#                               top_content[index_of_forcefield + len('; Include forcefield parameters'):])
#    else:
#        updated_top_content = top_content + '\n' + itp_files_section
#    with open(top_file_path, 'w', encoding='utf-8') as top_file:
#        top_file.write(updated_top_content)

    molecule_data = {
        'Protein_A': 1,
        'SOL': 10832
    }
    md.update_top_file('./' + top_file_path, molecule_data)

    commands = [
        f'gmx grompp -f {minimization_mdp} -c {input_molecule}_solv.gro -p topol.top -o em.tpr -maxwarn 3',    # 에너지 최적화 준비
        f'gmx mdrun -v -deffnm em',  # 에너지 최적화 실행
        f'gmx grompp -f {nvt_mdp} -c em.gro -p topol.top -o nvt.tpr -maxwarn 4', # NVT 준비
        f'gmx mdrun -v -deffnm nvt', # NVT 실행
        f'gmx grompp -f {npt_mdp} -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 4', # NPT 준비
        f'gmx mdrun -v -deffnm npt', # NPT 실행
        f'gmx grompp -f {mdp_file_path} -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr -maxwarn 3', # MD 준비
        f'gmx mdrun -v -deffnm md_0_1' # MD 실행
    ]

    for command in commands:
        md.run_gmx_command(command, "작업 완료")
        md.print_rotating_bar()

    # 사용자 입력 및 추가 계산 실행
    user_input = input("계산할 것을 선택하세요. (1: Energy, 2: RMSD, 3: RMSF, 4: Gyrate) : ")
    input_list = user_input.split()

    for option in input_list:
        if option == '1':
            command = f'gmx energy -f md_0_1.edr -o potential.xvg'
            md.run_gmx_command(command, "에너지 계산 완료")
            md.print_rotating_bar()
            md.print_rotating_bar()
            md.print_rotating_bar()

        elif option == '2':
            command = f'gmx rms -s md_0_1.tpr -f md_0_1.xtc -o rmsd.xvg'
            md.run_gmx_command(command, "RMSD 계산")
            md.print_rotating_bar()
            md.print_rotating_bar()
            md.print_rotating_bar()
        elif option == '3':
            command = f'gmx rmsf -s md_0_1.tpr -f md_0_1.xtc -o rmsf.xvg'
            md.run_gmx_command(command, "RMSF 계산 완료")
            md.print_rotating_bar()
            md.print_rotating_bar()
            md.print_rotating_bar()
        elif option == '4':
            command = f'gmx gyrate -s md_0_1.tpr -f md_0_1.xtc -o gyrate.xvg'
            md.run_gmx_command(command, "Gyrate 계산 완료")
            md.print_rotating_bar()
            md.print_rotating_bar()
            md.print_rotating_bar()

    md.cleanup_files()
    print("Job Success!")

if __name__ == "__main__":
    main()