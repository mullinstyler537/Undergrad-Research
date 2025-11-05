import cv2
import numpy as np

# -----------------------------
# Load SEM image
# -----------------------------
image_path = r"C:\Users\D00456326\Desktop\Research Photos\Nanoparticles\NOCOB2.JPG"
img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
if img is None:
    print("Error: Image not found.")
    exit()

# -----------------------------
# Ask for micron scale
# -----------------------------
micron_value = float(input("Enter the micron scale represented by the scale bar: "))
pixel_length = float(input("Enter the number of pixels corresponding to that scale: "))
micron_per_pixel = micron_value / pixel_length
print(f"Micron per pixel ratio: {micron_per_pixel:.6f} μm/pixel")

# -----------------------------
# Crop bottom 10%
# -----------------------------
height = img.shape[0]
img_cropped = img[:int(height * 0.9), :]

# -----------------------------
# Measuring Pass: preprocessing
# -----------------------------
# Minimal blur to preserve edges
img_blur = cv2.GaussianBlur(img_cropped, (3, 3), 0)

# Adaptive thresholding (invert so particles are white)
img_thresh = cv2.adaptiveThreshold(
    img_blur, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
    cv2.THRESH_BINARY_INV, 11, 2
)

# Morphological opening to remove tiny noise
kernel = np.ones((2, 2), np.uint8)
img_thresh = cv2.morphologyEx(img_thresh, cv2.MORPH_OPEN, kernel)

# -----------------------------
# Watershed to separate touching particles
# -----------------------------
dist_transform = cv2.distanceTransform(img_thresh, cv2.DIST_L2, 5)
_, sure_fg = cv2.threshold(dist_transform, 0.55 * dist_transform.max(), 255, 0)

sure_fg = np.uint8(sure_fg)
unknown = cv2.subtract(img_thresh, sure_fg)

_, markers = cv2.connectedComponents(sure_fg)
markers = markers + 1
markers[unknown == 255] = 0

img_color = cv2.cvtColor(img_cropped, cv2.COLOR_GRAY2BGR)
markers = cv2.watershed(img_color, markers)
img_thresh[markers == -1] = 0

# -----------------------------
# Find contours for measuring
# -----------------------------
contours, _ = cv2.findContours(img_thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
min_area_pixels = 5  # remove tiny noise
contours = [c for c in contours if cv2.contourArea(c) >= min_area_pixels]

areas_pixels = [cv2.contourArea(c) for c in contours]
areas_microns = [a * (micron_per_pixel ** 2) for a in areas_pixels]

# -----------------------------
# Output average area
# -----------------------------
if len(areas_microns) == 0:
    print("No nanoparticles detected for measuring.")
else:
    avg_area = np.mean(areas_microns)
    print(f"Average nanoparticle area (measuring pass): {avg_area:.4f} μm²")

# Optional: save visualization
cv2.imwrite("measuring_pass.png", img_thresh)
print("Thresholded image for measuring saved as 'measuring_pass.png'.")
