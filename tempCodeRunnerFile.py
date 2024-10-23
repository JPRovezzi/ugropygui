def show_molecule_svg(molecule_name):
    mol = Chem.MolFromSmiles("CCCCCC")
    if mol:
        svg = Draw.MolsToGridImage([mol], useSVG=True)
        return svg
    else:
        return None

svg = show_molecule_svg(molecule_name)
if svg:
    svg_window = tk.Tk()
    svg_label = tk.Label(svg_window, text="Molecule Structure")
    svg_label.pack()
    
    svg_image = tk.Label(svg_window, image=svg.data)
    svg_image.pack()
    
    close_button = tk.Button(svg_window, text="Close", command=svg_window.destroy)
    close_button.pack()
    center_window(svg_window)
    svg_window.mainloop()
else:
    print("Invalid molecule structure.")