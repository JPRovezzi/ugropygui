# Description: This script is the main script for the UGROpyGUI application. 
# It takes a molecule name as input and displays the molecule groups in a GUI window.

# Import the required libraries
# UgroPy is a Python library that provides a simple interface to the UNIFAC group contribution method.
from ugropy import Groups
from ugropy import unifac
# Tkinter is a standard GUI library for Python.
import tkinter as tk
from tkinter import simpledialog
# The subprocess library is used to run external commands.
import subprocess
# The os library provides a way to interact with the operating system.
import os
# The PIL library is used to work with images.
from PIL import Image, ImageDraw, ImageFont

#------------------------------------------------------------
# Functions to display the GUI windows

def window_welcome():
    '''Display a welcome window.'''
    window = tk.Tk()
    window.title(title)
    label = tk.Label(window, text="\nWelcome to UGROpyGUI!\n")
    label.pack()
    start_button = tk.Button(window, text="Start", command=window.destroy)
    start_button.pack()
    
    set_window_size(window)
    center_window(window)
    window.mainloop()

def window_selection():
    '''Display a window to select the input type (name or SMILES code).'''
    window = tk.Tk()
    window.title(title)
    label = tk.Label(window, text="\nSelect the input type:\n")
    label.pack()
    global input_type
    input_type = None
    def select_name():
        window.destroy()
        global input_type
        input_type = 'name'
    
    def select_smiles():
        window.destroy()
        global input_type
        input_type = 'smiles'

    def select_exit():
        window.destroy()
        global input_type
        input_type = 'exit'
    
    name_button = tk.Button(window, text="Name", command=select_name)
    name_button.pack()
    
    smiles_button = tk.Button(window, text="SMILES", command=select_smiles)
    smiles_button.pack()
    label = tk.Label(window, text="\n")
    label.pack()
    smiles_button = tk.Button(window, text="Exit", command=select_exit)
    smiles_button.pack()
    
    
    set_window_size(window)
    center_window(window)
    window.mainloop()
    print(input_type)
    return input_type

def window_get_name_old(text):
    '''Display a window to get the molecule name from the user.'''
    root = tk.Tk()
    root.title(title)
    root.withdraw()  # Hide the main window
    molecule_name = simpledialog.askstring("Input", "Please enter the molecule "+text+" :")
    
    set_window_size(root)
    center_window(root)
    root.destroy()
    return molecule_name
def window_get_name(text):
    '''Display a window to get the molecule name from the user.'''
    window = tk.Tk()
    window.title(title)
    label = tk.Label(window, text="Please enter the molecule " + text + " :")
    label.pack(pady=10)
    molecule_name_var = tk.StringVar()
    entry = tk.Entry(window, textvariable=molecule_name_var)
    entry.pack(pady=10)
    def submit():
        window.destroy()
    submit_button = tk.Button(window, text="Submit", command=submit)
    submit_button.pack(pady=10)
    
    set_window_size(window)
    center_window(window)
    window.mainloop()
    return molecule_name_var.get()

def window_message(message):
    '''Display a window with a message.'''
    window = tk.Tk()
    window.title(title)
    label = tk.Label(window, text="\n"+message+"\n")
    #image = tk.PhotoImage(file="./output.png")
    #image_label = tk.Label(window, image=image)
    #image_label.pack()
    label.pack()
    close_button = tk.Button(window, text="Close", command=window.destroy)
    close_button.pack()
    
    set_window_size(window)
    center_window(window)
    window.mainloop()

def window_picture(message):
    '''Display a window with a message and an image.'''
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
    
    set_window_size(window)
    center_window(window)
    window.mainloop()

def center_window(window):
    '''Center the window on the screen.'''
    window.update_idletasks()
    width = window.winfo_width()
    height = window.winfo_height()
    x = (window.winfo_screenwidth() // 2) - (width // 2)
    y = (window.winfo_screenheight() // 2) - (height // 2)
    window.geometry(f'{width}x{height}+{x}+{y}')
    return None

def set_window_size(window, width=640, height=480):
    '''Set the window size to 640x480 pixels by default.'''
    window.geometry(f'{width}x{height}')
    return None

def run_c_script(c_script_path: str, output_path: str):
    '''Compile and run a C script.'''
    if not os.path.exists("SvgToPng.exe"):  
        # Compile the C script
        compile_command = f"gcc {c_script_path} -o {output_path}"
        compile_process = subprocess.run(compile_command, shell=True, check=True)
    # Run the compiled C script
    run_command = f"./{output_path}"
    #run_process = subprocess.run(run_command, shell=True, check=True, capture_output=True, text=True)
    run_process = subprocess.run(run_command)
    return run_process.stdout

def write_groups2picture(dict):
    '''Write the molecule groups in the picture.'''
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
title = "UGROpyGUI" # Set the title of the GUI window


window_welcome() # Display the welcome window    
while True:
    input_type = window_selection()
    if input_type != 'exit' and input_type != None:
        molecule_name = window_get_name(input_type) # Get the molecule name from the user
        if molecule_name:
            try:
                molecule = Groups(identifier=molecule_name,identifier_type=input_type)
                molecule_groups = unifac.get_groups(identifier=molecule_name,identifier_type=input_type)
            except:
                window_message("Invalid molecule name.")
                continue
            svg_string = molecule_groups.get_solution_svg() # Get the SVG information
            with open("input.svg", "w") as file:
                file.write(svg_string) #Save the SVG
            run_c_script("SvgToPng.c", "SvgToPng") # Run the C script to convert the SVG to PNG
            write_groups2picture(molecule.unifac.subgroups) # Write the molecule groups in the picture
            window_picture("The molecule groups are displayed below.")

        else:
            print("No molecule name provided.")
            continue # Exit the loop if no molecule name is provided
    else:
        break # Exit the loop if no input type is selected



