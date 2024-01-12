import subprocess
import time
import sys
from datetime import datetime

def print_save_timestamp():
    save_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"Code saved at: {save_timestamp}")


def read_lines_from_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.readlines()

def extract_itp_files_names(itp_files):
    itp_files_name = []
    for itp_file in itp_files:
        itp_files_name.append(itp_file.split('/')[-1].strip())
    return itp_files_name


def run_gmx_command(command, message):
    print(f"Command : {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print(message)
        print ("====================================")
    except subprocess.CalledProcessError as e:
        print(f"오류 발생: {e}")
        print ("====================================")

def sleep(seconds):
    time.sleep(seconds)

def cleanup_files():
    subprocess.run('rm *.top', shell=True, check=False)
    subprocess.run('rm mdout.mdp', shell=True, check=False)
    subprocess.run('rm posre.itp', shell=True, check=False)
    subprocess.run('rm *.gro', shell=True, check=False)
    subprocess.run('rm *.top', shell=True, check=False)
    subprocess.run('rm #topol.top.2#', shell=True, check=False)
    subprocess.run('rm #topol.top.1#', shell=True, check=False)
    subprocess.run('rm *.trr', shell=True, check=False)
    subprocess.run('rm *.tpr', shell=True, check=False)
    subprocess.run('rm *.log', shell=True, check=False)
    

def print_rotating_bar(duration_seconds=3):
    rotation_symbols = ['\\', '|', '/', '-']
    start_time = time.time()

    while time.time() - start_time < duration_seconds:
        for symbol in rotation_symbols:
            sys.stdout.write('\r' + symbol)
            sys.stdout.flush()
            time.sleep(0.25)

    sys.stdout.write('\r')  # Clear the line
    sys.stdout.flush()



def update_top_file(top_file_path, molecule_data):
    try:
        # top 파일 읽기
        with open(top_file_path, 'r') as f:
            lines = f.readlines()

        # [ molecules ] 섹션 찾기
        molecules_start = lines.index('[ molecules ]\n')
        molecules_end = lines.index('\n', molecules_start)

        # [ molecules ] 섹션 업데이트
        lines[molecules_start + 1:molecules_end] = [f'{molecule} {count}\n' for molecule, count in molecule_data.items()]

        # top 파일 쓰기
        with open(top_file_path, 'w') as f:
            f.writelines(lines)

        print("topol.top 파일 업데이트 완료")
    except Exception as e:
        print(f"오류 발생: {e}")


def check_gmx_file(avail_file):
    try 


def additional_job():
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