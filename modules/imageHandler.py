
# Importing required libraries

# Tkinter is a standard GUI library for Python.
import tkinter as tk

# The PIL library is used to work with images.
from PIL import ImageTk

def insert_image(where, image_path):
    '''Insert an image into a frame.'''
    image = ImageTk.PhotoImage(file=image_path)
    image_widget = tk.Label(
        where,
        image = image
        )
    image_widget.image = image
    image_widget.pack()
    return None