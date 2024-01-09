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