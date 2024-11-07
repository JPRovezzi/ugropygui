# Description: This script is the main script for the UGROpyGUI application. 
# It takes a molecule name as input and displays the molecule groups in a GUI window.

# Import the required libraries
# UgroPy is a Python library that provides a simple interface to the UNIFAC group contribution method.
from ugropy import Groups
from ugropy import unifac
# Tkinter is a standard GUI library for Python.
import tkinter as tk

#from tkinter import simpledialog
# The subprocess library is used to run external commands.
import subprocess
# The os library provides a way to interact with the operating system.
import os
# The PIL library is used to work with images.
from PIL import Image, ImageDraw, ImageFont, ImageTk

#------------------------------------------------------------
# Functions to display the GUI windows

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
    #print("picture loaded")
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

def load_frame_welcome():
    clear_widgets_except(frame_welcome)
    frame_welcome.tkraise()
    frame_welcome.pack_propagate(False)
    # frame_welcome widgets
    tk.Label(
        frame_welcome, 
        text="\nWelcome to UGROpyGUI!\n",
        bg=bg_color,
        fg="white",
        font=("TkMenuFont",14)
        ).pack()

    tk.Button(
        frame_welcome,
        text="START",
        fg="black",
        font=("TkMenuFont",12),
        bg="white",
        cursor="hand2",
        activebackground="gray",
        command=lambda:load_frame_selection()
        ).pack(pady=10)
    
    tk.Button(
    frame_welcome,
    text="EXIT",
    fg="black",
    font=("TkMenuFont",12),
    bg="white",
    cursor="hand2",
    activebackground="gray",
    command=lambda:root.destroy()
    ).pack(pady=10)

    return None

def load_frame_selection(error_message=None):
    def select_name():
        #global input_type
        input_type = 'name'
        molecule_id = molecule_id_var.get()
        outcome = get_results(molecule_id,input_type)
        molecule,error = outcome
        if error == None and molecule != None:
            load_frame_result(molecule)
            error_message = ""
        elif error == 1:
            clear_widgets_except(None)
            error_message = "The NAME identifier is not valid. Please try again."
            load_frame_selection(error_message)
        elif error == 2:
            clear_widgets_except(None)
            error_message = "You must enter a NAME identifier. Please try again."
            load_frame_selection(error_message)
        else: 
            error_message = "Unknown error. Please try again."
            clear_widgets_except(None)
            load_frame_selection(error_message)

    def select_smiles():
        #global input_type
        input_type = 'smiles'
        
        molecule_id = molecule_id_var.get()
        outcome = get_results(molecule_id,input_type)
        molecule,error = outcome
        if error == None and molecule != None:
            load_frame_result(molecule)
            error_message = ""
        elif error == 1:
            clear_widgets_except(None)
            error_message = "The SMILES identifier is not valid. Please try again."
            load_frame_selection(error_message)
        elif error == 2:
            clear_widgets_except(None)
            error_message = "You must enter a SMILES identifier. Please try again."
            load_frame_selection(error_message)
        else: 
            error_message = "Unknown error. Please try again."
            clear_widgets_except(None)
            load_frame_selection(error_message)

    clear_widgets_except(frame_selection)
    frame_selection.tkraise()

    # frame_selection widgets
    label = tk.Label(frame_selection,text="",bg=bg_color).pack(pady=0)
    label = tk.Label(
        frame_selection, 
        text = "Please enter the Chemical identifier of the molecule:",
        bg=bg_color,
        fg="white",
        font=("TkMenuFont",14)
        ).pack(pady=0)
    molecule_id_var = tk.StringVar()
    entry = tk.Entry(
        frame_selection, 
        textvariable = molecule_id_var
        ).pack(pady=20)

    label = tk.Label(
        frame_selection, 
        text = error_message,
        bg=bg_color,
        fg="red",
        font=(14)
        ).pack(pady=0)

    tk.Label(
        frame_selection, 
        text="Select the type of Chemical identifier",
        bg=bg_color,
        fg="white",
        font=("TkMenuFont",14)
        ).pack()

    tk.Button(
        frame_selection,
        text="NAME",
        fg="black",
        font=("TkMenuFont",12),
        bg="white",
        cursor="hand2",
        activebackground="gray",
        command=lambda:select_name()
        ).pack(pady=10)

    tk.Button(
        frame_selection,
        text="SMILES",
        fg="black",
        font=("TkMenuFont",12),
        bg="white",
        cursor="hand2",
        activebackground="gray",
        command=lambda:select_smiles()
        ).pack(pady=10)

    tk.Button(
        frame_selection,
        text="BACK",
        fg="black",
        font=("TkMenuFont",12),
        bg="white",
        cursor="hand2",
        activebackground="gray",
        command=lambda:load_frame_welcome()
        ).pack(pady=50)

    return None

def load_frame_result(molecule):
    
    clear_widgets_except(frame_result)
    frame_result.tkraise()
    
    tk.Label(
        frame_result,
        text="The molecule groups are displayed below.",
        bg=bg_color,
        fg="white",
        font=("TkMenuFont",14)
        ).pack(pady=20)
    
    image=ImageTk.PhotoImage(file="./output.png")
    image_widget = tk.Label(
        frame_result,
        image = image
        )
    image_widget.image = image
    image_widget.pack()

    tk.Label(
        frame_result,
        text = molecule.unifac.subgroups,
        bg=bg_color,
        fg="white",
        ).pack()

    tk.Button(
        frame_result,
        text="BACK",
        fg="black",
        font=("TkMenuFont",12),
        bg="white",
        cursor="hand2",
        activebackground="gray",
        command=lambda:load_frame_welcome()
        ).pack(pady=50)
    
    return None


def clear_widgets_except(currentFrame):
    for frame in frames:
        if frame != currentFrame:
            for widget in frame.winfo_children():
                widget.destroy()
    return None

def get_results(molecule_id,input_type):  
    outcome = (None,None)
    if molecule_id != None and molecule_id != "":
            try:
                #print(input_type, molecule_id)
                molecule = Groups(
                    identifier = molecule_id,
                    identifier_type = input_type
                    )
                molecule_groups = unifac.get_groups(
                    identifier = molecule_id,
                    identifier_type = input_type
                    )
            except:
                error = 1
                #print(1)
                outcome = (None, error)
                return outcome
            svg_string = molecule_groups.get_solution_svg() # Get the SVG information
            with open("input.svg", "w") as file:
                file.write(svg_string) #Save the SVG
            run_c_script("SvgToPng.c", "SvgToPng") # Run the C script to convert the SVG to PNG
            write_groups2picture(molecule.unifac.subgroups) # Write the molecule groups in the picture
            #print("The molecule groups are displayed below.")
            error = None
            #print(0)
            outcome = (molecule, error)
            return outcome

    else:
        error = 2
        #print(2)
        outcome = (None,error)
        return outcome
    
#------------------------------------------------------------

# Basic configuration
bg_color = '#f36c31' #IPQA Color

# Global variables
molecule_id = "" # Variable to store the molecule identifier

# Create the main window
root = tk.Tk()
root.title("UGROpyGUI")
root.resizable(0,0) # Disable resizing
root.eval("tk::PlaceWindow . center")

# Create a frame widget
frame_welcome = tk.Frame(root, width=640, height=480, bg=bg_color)
frame_selection = tk.Frame(root, width=640, height=480, bg=bg_color)
frame_getName = tk.Frame(root, width=640, height=480, bg=bg_color)
frame_result = tk.Frame(root, width=640, height=480, bg=bg_color)
frames = (frame_welcome,frame_selection, frame_getName, frame_result)

for frame in frames:
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    frame.grid(row=0, column=0, sticky="nesw")

load_frame_welcome()

root.mainloop()





