from ugropy import Groups
from ugropy import unifac
import tkinter as tk
from tkinter import simpledialog
from rdkit import Chem
from ugropy import DefaultSolver
from rdkit.Chem import Draw
from IPython.display import SVG
import subprocess
import os
from PIL import Image, ImageDraw, ImageFont

#------------------------------------------------------------

def get_molecule_name():
    root = tk.Tk()
    root.title(title)
    root.withdraw()  # Hide the main window
    molecule_name = simpledialog.askstring("Input", "Please enter the molecule name:")
    center_window(root)
    root.destroy()
    return molecule_name

def window_message(message):
    window = tk.Tk()
    window.title(title)
    label = tk.Label(window, text="\n"+message+"\n")
    #image = tk.PhotoImage(file="./output.png")
    #image_label = tk.Label(window, image=image)
    #image_label.pack()
    label.pack()
    close_button = tk.Button(window, text="Close", command=window.destroy)
    close_button.pack()
    center_window(window)
    window.mainloop()

def window_picture(message):
    window = tk.Tk()
    window.title(title)
    label = tk.Label(window, text="\n"+message+"\n")
    label.pack()
    label = tk.Label(window, text=molecule.unifac.subgroups)
    image = tk.PhotoImage(file="./output.png")
    image_label = tk.Label(window, image=image)
    image_label.pack()
    label.pack()
    close_button = tk.Button(window, text="Close", command=window.destroy)
    close_button.pack()
    center_window(window)
    window.mainloop()

def center_window(window):
    window.update_idletasks()
    width = window.winfo_width()
    height = window.winfo_height()
    x = (window.winfo_screenwidth() // 2) - (width // 2)
    y = (window.winfo_screenheight() // 2) - (height // 2)
    window.geometry(f'{width}x{height}+{x}+{y}')
    return None

def run_c_script(c_script_path: str, output_path: str):
    if not os.path.exists("SvgToPng.exe"):  
        # Compile the C script
        compile_command = f"gcc {c_script_path} -o {output_path}"
        compile_process = subprocess.run(compile_command, shell=True, check=True)
    # Run the compiled C script
    run_command = f"./{output_path}"
    #run_process = subprocess.run(run_command, shell=True, check=True, capture_output=True, text=True)
    run_process = subprocess.run(run_command)
    return run_process.stdout

def writeInPicture(dict):
    # Load the image
    image = Image.open('output.png')
    print("picture loaded")
    # Initialize ImageDraw
    draw = ImageDraw.Draw(image)

    # Define the text and position
    #text = list(dict.keys())[0]
    position = (25, 10)

    # Load a font
    font = ImageFont.load_default()

    # Add text to image
    for i, key in enumerate(dict.keys()):
        position = (25, 10 + i * 25)  # Adjust the position for each key
        draw.text(position, key, font=font, fill=(0, 0, 0))
    #draw.text(position, text, font=font, fill=(0, 0, 0))

    # Save the image
    image.save('output.png')
    return None

#------------------------------------------------------------
title = "UGROpyGUI"

while True:
    molecule_name = get_molecule_name()

    if molecule_name:
        try:
            molecule = Groups(molecule_name)
            molecule_groups = unifac.get_groups(molecule_name)
        except:
            print("Invalid molecule name.")
            window_message("Invalid molecule name.")
            continue
        svg_string = molecule_groups.get_solution_svg()
        with open("input.svg", "w") as file:
            file.write(svg_string)
        
        run_c_script("SvgToPng.c", "SvgToPng")
        writeInPicture(molecule.unifac.subgroups)
        window_picture("The molecule groups are displayed below.")
        print(molecule.unifac.subgroups)


    else:
        print("No molecule name provided.")
        exit()



