import tkinter as tk
import modules.guiFrames.functions as functions
import modules.guiFrames.frameRoot as frameRoot

def load():
        save_window = tk.Toplevel(frameRoot.root)
        save_window.resizable(0, 0)  # Disable resizing
        save_window.title("Save Last Picture")
        save_window.transient(frameRoot.root)
        save_window.grab_set()
        #save_window.geometry("400x300")
        
        # Display the output PNG as a widget
        
        img = tk.PhotoImage(file="output.png")
        img_label = tk.Label(save_window, image=img)
        img_label.image = img  # Keep a reference to avoid garbage collection
        img_label.pack(pady=10)
        
        # Entry box to write the path to save
        #tk.Label(save_window, text="Save Path:").pack(pady=5)
        #path_entry_text = tk.StringVar()
        #path_entry = tk.Entry(save_window, width=50, textvariable=path_entry_text)
        #path_entry_text.set("C:/")
        #path_entry.pack(pady=5)
        
        # Selection box to select if it is a SVG image or a PNG image
        tk.Label(save_window, text="Select Format:").pack(pady=5)
        format_var = tk.StringVar(value="png")
        
        # Frame to hold the radio buttons
        radioButton_frame = tk.Frame(save_window)
        radioButton_frame.pack(pady=10)

        tk.Radiobutton(
                radioButton_frame, 
                text="PNG", 
                variable = format_var, 
                value="png"
                ).pack(side=tk.LEFT,padx=0)
        tk.Radiobutton(
                radioButton_frame, 
                text="SVG", 
                variable = format_var, 
                value = "svg"
                ).pack(side=tk.LEFT,padx=0)
               
        # Frame to hold the buttons
        button_frame = tk.Frame(save_window)
        button_frame.pack(pady=10)

        # Save button
        tk.Button(
            button_frame, 
            text="Save", 
            command=lambda: functions.save_image(format_var.get())
            ).pack(side=tk.LEFT, padx=5)

        # Cancel button
        tk.Button(
            button_frame, 
            text="Close", 
            command=lambda: save_window.destroy()
            ).pack(side=tk.LEFT, padx=5)
