import subprocess
import time
import md_automation as md

def main():
    file_path = './input_file_path.txt'
    lines = md.read_lines_from_file(file_path)

    top_file_path = "./topol.top"
    pdb_file_path = lines[6].strip()
    mdp_file_path = lines[9].strip('\n')
    itp_files = lines[12:]

    itp_files_name = md.extract_itp_files_names(itp_files)

    md.print_header()
    print(f"Load pdb file : {pdb_file_path}")
    print(f"Load mdp file : {mdp_file_path}")
    print("====================================")

    input_molecule = pdb_file_path.replace('.pdb', '')
    print(f"input_molecule : {input_molecule}")

    # pdb2gmx 실행
    pdb2gmx_command = f'gmx pdb2gmx -f {input_molecule}.pdb -o {input_molecule}_processed.gro -water spce -ff oplsaa -p topol.top'
    md.run_gmx_command(pdb2gmx_command, "pdb2gmx 실행 완료")

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

    # grompp 및 MD 시뮬레이션 실행
    commands = [
        f'gmx editconf -f {input_molecule}_processed.gro -o {input_molecule}_newbox.gro -c -d 1.0 -bt cubic',
        f'gmx solvate -cp {input_molecule}_newbox.gro -cs spc216.gro -o {input_molecule}_solv.gro -p topol.top',
        f'gmx grompp -f {mdp_file_path} -c {input_molecule}_solv.gro -p topol.top -o em.tpr',
        f'gmx mdrun -v -deffnm em',
        f'gmx grompp -f {mdp_file_path} -c em.gro -p topol.top -o nvt.tpr',
        f'gmx mdrun -v -deffnm nvt',
        f'gmx grompp -f {mdp_file_path} -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr',
        f'gmx mdrun -v -deffnm npt',
        f'gmx grompp -f {mdp_file_path} -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr',
        f'gmx mdrun -v -deffnm md_0_1'
    ]

    for command in commands:
        md.run_gmx_command(command, "작업 완료")
        time.sleep(1)

    # 사용자 입력 및 추가 계산 실행
    user_input = input("계산할 것을 선택하세요. (1: 에너지, 2: RMSD, 3: RMSF, 4: Gyrate) : ")
    input_list = user_input.split()

    for option in input_list:
        if option == '1':
            command = f'gmx energy -f md_0_1.edr -o potential.xvg'
            md.run_gmx_command(command, "에너지 계산 완료")
            time.sleep(10)
        elif option == '2':
            command = f'gmx rms -s md_0_1.tpr -f md_0_1.xtc -o rmsd.xvg'
            md.run_gmx_command(command, "RMSD 계산")
            time.sleep(10)
        elif option == '3':
            command = f'gmx rmsf -s md_0_1.tpr -f md_0_1.xtc -o rmsf.xvg'
            md.run_gmx_command(command, "RMSF 계산 완료")
            time.sleep(10)
        elif option == '4':
            command = f'gmx gyrate -s md_0_1.tpr -f md_0_1.xtc -o gyrate.xvg'
            md.run_gmx_command(command, "Gyrate 계산 완료")
            time.sleep(10)

    md.cleanup_files()
    print("Job Success!")

if __name__ == "__main__":
    main()