from tkinter import Tk, filedialog
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- File picker ---
Tk().withdraw()
image_path = filedialog.askopenfilename(
    title="Select SEM image",
    filetypes=[("Image files", "*.jpg *.jpeg *.png *.tif *.tiff")]
)

img = cv2.imread(image_path)
if img is None:
    raise FileNotFoundError(f"Could not read image at: {image_path}")

# --- Ask user for scale information ---
micron_scale = float(input("Enter the scale bar length in microns: "))
pixel_scale = img.shape[0]

# Convert to nanometers per pixel
nm_per_pixel = (micron_scale * 1000) / pixel_scale

# --- Preprocessing ---
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
blur = cv2.GaussianBlur(gray, (5, 5), 0)
thresh = cv2.adaptiveThreshold(
    blur, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
    cv2.THRESH_BINARY_INV, 11, 2
)

# --- Contour detection ---
contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
areas_px = []

for cnt in contours:
    area = cv2.contourArea(cnt)
    if area > 10:
        areas_px.append(area)
        cv2.drawContours(img, [cnt], -1, (0, 0, 255), 1)

# --- Area calculations ---
areas_nm2 = np.array(areas_px) * (nm_per_pixel ** 2)
areas_um2 = areas_nm2 / 1e6

# --- Statistics ---
avg_area = np.mean(areas_um2)
std_area = np.std(areas_um2)

# --- Outputs ---
cv2.imwrite("sem_outlined.jpg", img)
pd.DataFrame(areas_um2, columns=["Area_µm²"]).to_csv("nanoparticle_areas.csv", index=False)

print(f"\nDetected {len(areas_um2)} nanoparticles")
print(f"Average area: {avg_area:.4f} µm²")
print(f"Std. dev.:    {std_area:.4f} µm²")

# --- Histogram ---
plt.figure(figsize=(6, 4))
plt.hist(areas_um2, bins=30, color='steelblue', edgecolor='black')
plt.xlabel("Particle area (µm²)")
plt.ylabel("Count")
plt.title("Nanoparticle Size Distribution")
plt.tight_layout()
plt.savefig("nanoparticle_size_distribution.png", dpi=300)
plt.show()
