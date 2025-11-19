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
micron_scale = float(input("Enter the vertical field-of-view (microns) for the full image height: "))

# Image height in pixels
pixel_scale = img.shape[0]

# Convert to nm per pixel
nm_per_pixel = (micron_scale * 1000) / pixel_scale
print(f"\nImage height: {pixel_scale} px")
print(f"Scale: {nm_per_pixel:.3f} nm/pixel")

# --- Preprocessing & segmentation ---
gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
blur = cv2.GaussianBlur(gray, (3, 3), 0)

thresh = cv2.adaptiveThreshold(
    blur, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
    cv2.THRESH_BINARY_INV, 11, 2
)

kernel = np.ones((3, 3), np.uint8)
opening = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel, iterations=1)
closing = cv2.morphologyEx(opening, cv2.MORPH_CLOSE, kernel, iterations=1)

sure_bg = cv2.dilate(closing, kernel, iterations=2)
dist = cv2.distanceTransform(closing, cv2.DIST_L2, 5)
_, sure_fg = cv2.threshold(dist, 0.35 * dist.max(), 255, 0)
sure_fg = np.uint8(sure_fg)
unknown = cv2.subtract(sure_bg, sure_fg)

_, markers = cv2.connectedComponents(sure_fg)
markers = markers + 1
markers[unknown == 255] = 0
markers = cv2.watershed(img, markers)

# --- Collect contours and calculate areas ---
areas_px = []
centroids = []
height, width = img.shape[:2]

for label in np.unique(markers):
    if label <= 1:
        continue
    mask = np.uint8(markers == label) * 255
    cnts, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    for cnt in cnts:
        area = cv2.contourArea(cnt)
        if area > 5:  # ignore very small noise
            # Compute centroid
            M = cv2.moments(cnt)
            if M["m00"] == 0:
                continue
            cX = int(M["m10"] / M["m00"])
            cY = int(M["m01"] / M["m00"])
            
            # Ignore particles touching image edges
            margin = 2  # pixels
            if cX < margin or cX > width - margin or cY < margin or cY > height - margin:
                continue
            
            areas_px.append(area)
            centroids.append((cX, cY))
            cv2.drawContours(img, [cnt], -1, (0, 0, 255), 1)

areas_px = np.array(areas_px)

# --- Automatic filtering based on median ± factor ---
median_area = np.median(areas_px)
lower_limit = 0.2 * median_area
upper_limit = 3 * median_area  # merged particles often >3× median
filtered_areas_px = areas_px[(areas_px >= lower_limit) & (areas_px <= upper_limit)]

# --- Area calculations ---
areas_nm2 = filtered_areas_px * (nm_per_pixel ** 2)
areas_um2 = areas_nm2 / 1e6
areas_diameter_nm = 2 * np.sqrt(areas_nm2 / np.pi)

# --- Statistics ---
avg_area = np.mean(areas_um2)
std_area = np.std(areas_um2)
avg_diameter = np.mean(areas_diameter_nm)
std_diameter = np.std(areas_diameter_nm)

# --- Outputs ---
cv2.imwrite("sem_outlined_filtered.jpg", img)
df = pd.DataFrame({
    "Area_µm²": areas_um2,
    "Diameter_nm": areas_diameter_nm
})
df.to_csv("nanoparticle_areas_filtered.csv", index=False)

print(f"\nDetected {len(areas_um2)} nanoparticles after filtering")
print(f"Average area: {avg_area:.4f} µm² (std: {std_area:.4f})")
print(f"Average diameter: {avg_diameter:.2f} nm (std: {std_diameter:.2f})")

# --- Histogram ---
plt.figure(figsize=(6, 4))
plt.hist(areas_um2, bins=30, color='steelblue', edgecolor='black')
plt.xlabel("Particle area (µm²)")
plt.ylabel("Count")
plt.title("Nanoparticle Size Distribution (Filtered)")
plt.tight_layout()
plt.savefig("nanoparticle_size_distribution_filtered.png", dpi=300)
plt.show()
