from PIL import Image, ImageDraw, ImageFont

# Load the image
image = Image.open('path_to_your_image.png')

# Initialize ImageDraw
draw = ImageDraw.Draw(image)

# Define the text and position
text = "Hola"
position = (20, 20)

# Load a font
font = ImageFont.load_default()

# Add text to image
draw.text(position, text, font=font, fill=(255, 255, 255))

# Save the image
image.save('path_to_save_image.png')