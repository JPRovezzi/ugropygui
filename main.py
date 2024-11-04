from ugropy import Groups
from ugropy import unifac
import tkinter as tk
from tkinter import simpledialog
from rdkit import Chem
from ugropy import DefaultSolver
from rdkit.Chem import Draw
from IPython.display import SVG
import subprocess

#------------------------------------------------------------

def get_molecule_name():
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    molecule_name = simpledialog.askstring("Input", "Please enter the molecule name:")
    center_window(root)
    root.destroy()
    return molecule_name

def center_window(window):
    window.update_idletasks()
    width = window.winfo_width()
    height = window.winfo_height()
    x = (window.winfo_screenwidth() // 2) - (width // 2)
    y = (window.winfo_screenheight() // 2) - (height // 2)
    window.geometry(f'{width}x{height}+{x}+{y}')

def run_c_script(c_script_path: str, output_path: str):
    # Compile the C script
    compile_command = f"gcc {c_script_path} -o {output_path}"
    compile_process = subprocess.run(compile_command, shell=True, check=True)
    
    # Run the compiled C script
    run_command = f"./{output_path}"
    #run_process = subprocess.run(run_command, shell=True, check=True, capture_output=True, text=True)
    run_process = subprocess.run(run_command)
    return run_process.stdout

#------------------------------------------------------------

molecule_name = get_molecule_name()

    
if molecule_name:
    try:
        molecule_groups = unifac.get_groups(molecule_name)
    except:
        print("Invalid molecule name.")
        exit()
    svg_string = molecule_groups.get_solution_svg()
    with open("input.svg", "w") as file:
        file.write(svg_string)
    molecule = Groups(molecule_name)
    c_script_path = "SvgToPng.c"
    output_path = "SvgToPng"
    output = run_c_script(c_script_path, output_path)

    goodbye_window = tk.Tk()
    goodbye_label = tk.Label(goodbye_window, text=molecule.unifac.subgroups)
    image = tk.PhotoImage(file="./output.png")
    image_label = tk.Label(goodbye_window, image=image)
    image_label.pack()
    goodbye_label.pack()
    close_button = tk.Button(goodbye_window, text="Close", command=goodbye_window.destroy)
    close_button.pack()
    center_window(goodbye_window)
    goodbye_window.mainloop()
    print(molecule.unifac.subgroups)


else:
    print("No molecule name provided.")




