
"""
    This module provides functions to create and manage a graphical user interface (GUI) for UGROpyGUI.
    Functions:
        load_frame_welcome(): Load the welcome frame with its widgets.
"""
# Import the required libraries:

# Tkinter is a standard GUI library for Python.
import tkinter as tk

# widgetClasses is a module that provides classes for GUI widgets.
import modules.widgetClasses as widgetClasses

#
import modules.guiFrames.frameRoot as frameRoot
import modules.guiFrames.frameSelection as frameSelection


#------------------------------------------------------------
import modules.guiFrames.functions as functions

def load():
    functions.clear_widgets_except(frameRoot.frameWelcome,frameRoot.frames)
    frameRoot.frame_welcome.tkraise()
    frameRoot.frame_welcome.pack_propagate(False)
    # frame_welcome widgets
    
    widgetClasses.TitleLabel(
        frameRoot.frame_welcome,
        text="\nWelcome to UGROpyGUI!\n"
        ).pack(pady=0)

    tk.Button(
        frameRoot.frame_welcome,
        text="START",
        fg="black",
        font=("TkMenuFont",12),
        bg="white",
        cursor="hand2",
        activebackground="gray",
        command=lambda:frameSelection.load()
        ).pack(pady=10)
    
    tk.Button(
    frameRoot.frame_welcome,
    text="EXIT",
    fg="black",
    font=("TkMenuFont",12),
    bg="white",
    cursor="hand2",
    activebackground="gray",
    command=lambda:frameRoot.root.destroy()
    ).pack(pady=10)

    return None