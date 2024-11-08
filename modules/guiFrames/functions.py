def clear_widgets_except(currentFrame,frames):
    for frame in frames:
        if frame != currentFrame:
            for widget in frame.winfo_children():
                widget.destroy()
    return None
