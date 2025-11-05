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
# Crop bottom 10%
# -----------------------------
height = img.shape[0]
img_cropped = img[:int(height * 0.9), :]

# -----------------------------
# Counting Pass: preprocessing
# -----------------------------
# Stronger blur to merge tiny fragments
img_blur = cv2.GaussianBlur(img_cropped, (7, 7), 0)

# Adaptive thresholding (invert so particles are white)
img_thresh = cv2.adaptiveThreshold(
    img_blur, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
    cv2.THRESH_BINARY_INV, 11, 2
)

# Morphological closing to merge touching particles
kernel = np.ones((3, 3), np.uint8)
img_thresh = cv2.morphologyEx(img_thresh, cv2.MORPH_CLOSE, kernel)

# Morphological opening to remove tiny noise
img_thresh = cv2.morphologyEx(img_thresh, cv2.MORPH_OPEN, kernel)

# -----------------------------
# Find contours for counting
# -----------------------------
contours, _ = cv2.findContours(img_thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
# Remove very small noise fragments
min_area_pixels = 5
contours = [c for c in contours if cv2.contourArea(c) >= min_area_pixels]

total_count = len(contours)
print(f"Total nanoparticles detected (counting pass): {total_count}")

# Optional: save visualization
cv2.imwrite("counting_pass.png", img_thresh)
print("Thresholded image for counting saved as 'counting_pass.png'.")
