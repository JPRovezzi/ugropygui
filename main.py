from ugropy import Groups
import tkinter as tk
from tkinter import simpledialog
from rdkit import Chem
from ugropy import DefaultSolver
from rdkit.Chem import Draw
from IPython.display import SVG

#usermol = Groups("2-ethylhexyl acrylate")

#print(usermol.unifac.subgroups)
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

molecule_name = get_molecule_name()

if molecule_name:
    molecule = Groups(molecule_name)
    goodbye_window = tk.Tk()
    goodbye_label = tk.Label(goodbye_window, text=molecule.unifac.subgroups)

    goodbye_label.pack()
    #goodbye_window.after(2000, goodbye_window.destroy)  # Close after 2 seconds
    close_button = tk.Button(goodbye_window, text="Close", command=goodbye_window.destroy)
    close_button.pack()
    center_window(goodbye_window)
    goodbye_window.mainloop()
    print(molecule.unifac.subgroups)


else:
    print("No molecule name provided.")



