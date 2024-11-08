
"""
    This module provides functions to create and manage a graphical user interface (GUI) for UGROpyGUI.
    Functions:
        clear_widgets_except(currentFrame): Clear all widgets except those in the specified frame.
"""
def clear_widgets_except(currentFrame,frames):
    for frame in frames:
        if frame != currentFrame:
            for widget in frame.winfo_children():
                widget.destroy()
    return None
