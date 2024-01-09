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
    except subprocess.CalledProcessError as e:
        print(f"오류 발생: {e}")

def sleep(seconds):
    time.sleep(seconds)

def cleanup_files():
    subprocess.run('rm *.top', shell=True, check=False)
    subprocess.run('rm mdout.mdp', shell=True, check=False)
    subprocess.run('rm posre.itp', shell=True, check=False)
    subprocess.run('rm *.gro', shell=True, check=False)
    subprocess.run('rm *.top', shell=True, check=False)

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