
"""
    This module provides functions to create and manage a graphical user interface (GUI) for UGROpyGUI.
    Functions:
        create_gui(): Create the main GUI window and frames.
        load_frame_welcome(): Load the welcome frame with its widgets.
        load_frame_selection(error_message=None): Load the selection frame with its widgets.
        load_frame_result(molecule): Load the result frame with its widgets.
        clear_widgets_except(currentFrame): Clear all widgets except those in the specified frame.
"""
# Import the required libraries:

# Tkinter is a standard GUI library for Python.
import tkinter as tk
# Configparser is used to read configuration files.
import configparser
# SvgHandler is a module that provides functions to handle SVG files.
import modules.svgHandler as svgHandler
# The PIL library is used to work with images.
from PIL import ImageTk
# ImageHandler is a module that provides functions to handle images.
import modules.imageHandler as imageHandler
# widgetClasses is a module that provides classes for GUI widgets.
import modules.widgetClasses as widgetClasses
#
import modules.guiFrames.frameWelcome as frameWelcome

#------------------------------------------------------------
# Read configuration file
config = configparser.ConfigParser()
config.read('config.cfg')

# Get background color from configuration file
bg_color = config.get('Settings', 'bg_color')

#------------------------------------------------------------


#------------------------------------------------------------
# Function definitions

def create_gui():

    '''Create the main GUI window and frames.'''
    global root, frame_welcome, frame_selection, frame_getName, frame_result, frames
    root = tk.Tk()
    root.title("UGROpyGUI")
    root.resizable(0, 0)  # Disable resizing
    root.eval("tk::PlaceWindow . center")

    frame_welcome = tk.Frame(root, width=640, height=480, bg=bg_color)
    frame_selection = tk.Frame(root, width=640, height=480, bg=bg_color)
    frame_getName = tk.Frame(root, width=640, height=480, bg=bg_color)
    frame_result = tk.Frame(root, width=640, height=480, bg=bg_color)
    frames = (frame_welcome, frame_selection, frame_getName, frame_result)

    for frame in frames:
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        frame.grid(row=0, column=0, sticky="nesw")

    frameWelcome.load()
    root.mainloop()
    return None





    





    
    