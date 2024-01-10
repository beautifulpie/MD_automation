import tkinter as tk
from tkinter import filedialog, scrolledtext
import os
import subprocess
from tkinter import ttk

class MDGui(tk.Tk):
    def __init__(self):
        super().__init__()

        # Selected file paths
        self.selected_pdb_file = ""
        self.selected_md_mdp_file = ""
        self.selected_minim_mdp_file = ""
        self.selected_npt_mdp_file = ""
        self.selected_nvt_mdp_file = ""
        self.selected_itp_files = []

        # Progress value
        self.progress_value = 0

        # GUI setup
        self.title('MD Automation GUI')
        self.geometry('800x600')

        self.label = tk.Label(self, text='MD Automation', font=('Helvetica', 16, 'bold'))
        self.label.pack(pady=10)

        self.create_file_selection_button('Select PDB File', self.get_pdb_file)
        self.create_file_selection_button('Select md.mdp', self.get_md_mdp_file)
        self.create_file_selection_button('Select minim.mdp', self.get_minim_mdp_file)
        self.create_file_selection_button('Select npt.mdp', self.get_npt_mdp_file)
        self.create_file_selection_button('Select nvt.mdp', self.get_nvt_mdp_file)

        self.create_file_selection_button('Select ITP Files', self.get_itp_files)

        self.progress_bar = ttk.Progressbar(self, orient=tk.HORIZONTAL, length=100, mode='determinate')
        self.progress_bar.pack(pady=10)

        self.output_text = scrolledtext.ScrolledText(self, width=50, height=10, wrap=tk.WORD)
        self.output_text.pack(pady=10)

        self.run_button = tk.Button(self, text='Run MD Automation', command=self.run_md_automation, font=('Helvetica', 12))
        self.run_button.pack(pady=10)

    def create_file_selection_button(self, button_text, click_function):
        button = tk.Button(self, text=button_text, command=click_function, font=('Helvetica', 10))
        button.pack(pady=5)

    def get_pdb_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("PDB Files", "*.pdb"), ("All Files", "*.*")])

        if file_path:
            self.selected_pdb_file = file_path
            self.output_text.insert(tk.END, f"Selected PDB file: {self.selected_pdb_file}\n")

    def get_itp_files(self):
        file_paths = filedialog.askopenfilenames(filetypes=[("ITP Files", "*.itp"), ("All Files", "*.*")])

        if file_paths:
            self.selected_itp_files = file_paths
            self.output_text.insert(tk.END, f"Selected ITP files: {', '.join(self.selected_itp_files)}\n")

    def get_md_mdp_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("MDP Files", "*.mdp"), ("All Files", "*.*")])

        if file_path:
            self.selected_md_mdp_file = file_path
            self.output_text.insert(tk.END, f"Selected md.mdp file: {self.selected_md_mdp_file}\n")

    def get_minim_mdp_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("MDP Files", "*.mdp"), ("All Files", "*.*")])

        if file_path:
            self.selected_minim_mdp_file = file_path
            self.output_text.insert(tk.END, f"Selected minim.mdp file: {self.selected_minim_mdp_file}\n")

    def get_npt_mdp_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("MDP Files", "*.mdp"), ("All Files", "*.*")])

        if file_path:
            self.selected_npt_mdp_file = file_path
            self.output_text.insert(tk.END, f"Selected npt.mdp file: {self.selected_npt_mdp_file}\n")

    def get_nvt_mdp_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("MDP Files", "*.mdp"), ("All Files", "*.*")])

        if file_path:
            self.selected_nvt_mdp_file = file_path
            self.output_text.insert(tk.END, f"Selected nvt.mdp file: {self.selected_nvt_mdp_file}\n")

    def run_md_automation(self):
        inst = 0
        if not self.selected_pdb_file:
            self.output_text.insert(tk.END, "Please add pdb file.\n")
            inst = 1
        if not self.selected_md_mdp_file:
            self.output_text.insert(tk.END, "Please add md.mdp file.\n")
            inst = 1
        if not self.selected_minim_mdp_file:
            self.output_text.insert(tk.END, "Please add minim.mdp file.\n")
            inst = 1
        if not self.selected_npt_mdp_file:
            self.output_text.insert(tk.END, "Please add npt.mdp file.\n")
            inst = 1
        if not self.selected_nvt_mdp_file:
            self.output_text.insert(tk.END, "Please add nvt.mdp file.\n")
            inst = 1
        if inst == 1:
            return

        # Run MD automation with the selected files
        self.label.config(text="MD Automation is running... Please wait.")
        self.output_text.insert(tk.END, "====================================\n")
        self.output_text.insert(tk.END, "MD_automation\n")
        self.output_text.insert(tk.END, "Made by Jung Youngwoo\n")
        self.output_text.insert(tk.END, "====================================\n")
        self.output_text.insert(tk.END, "Generate the input_path_file\n")

        # Generating input path
        pdb_file_name = os.path.splitext(os.path.basename(self.selected_pdb_file))[0]
        minim_mdp_name = os.path.splitext(os.path.basename(self.selected_minim_mdp_file))[0]
        npt_mdp_name = os.path.splitext(os.path.basename(self.selected_npt_mdp_file))[0]
        nvt_mdp_name = os.path.splitext(os.path.basename(self.selected_nvt_mdp_file))[0]
        md_mdp_name = os.path.splitext(os.path.basename(self.selected_md_mdp_file))[0]
        itp_files_name = [os.path.splitext(os.path.basename(file))[0] for file in self.selected_itp_files]

        with open('임시.txt', 'w') as f:
            f.write("====================================\n")
            f.write("MD_automation\n")
            f.write("\n")
            f.write("Made by Youngwoo Jung\n")
            f.write("====================================\n")
            f.write(f"structure file name : \n {pdb_file_name}\n")
            f.write("\n")
            f.write(f"Simulator parameter file name : \n {md_mdp_name} \n {minim_mdp_name} \n {nvt_mdp_name} \n {npt_mdp_name}\n")
            f.write("\n")
            f.write("Itp file names : \n")
            for name in itp_files_name:
                f.write(f"{name}\n")

        print(f"pdb_file_path : {self.selected_pdb_file}")
        print(f'minim_mdp_path : {self.selected_minim_mdp_file}')
        print(f'npt_mdp_path : {self.selected_npt_mdp_file}')
        print(f'nvt_mdp_path : {self.selected_nvt_mdp_file}')
        print(f'md_mdp_path : {self.selected_md_mdp_file}')
        print(f'itp_files_path : {self.selected_itp_files}')
        print(f'pdb_file_name : {pdb_file_name}')
        print(f'minim_mdp_name : {minim_mdp_name}')

        self.progress_bar['maximum'] = 100
        self.progress_bar['value'] = 0

        # Running main.py
        subprocess.run('python main.py', shell=True, check=True)
        print("Run the main.py")
        print("\n\n")

        # Display success message
        self.label.config(text="MD Automation completed successfully!")
        self.output_text.insert(tk.END, "Job Success!\n")

if __name__ == '__main__':
    app = MDGui()
    app.mainloop()
