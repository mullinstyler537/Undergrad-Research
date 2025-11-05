import cv2
import numpy as np

def get_micron_to_pixel_ratio():
    """
    Prompt the user for a known distance in microns and the corresponding pixel distance,
    then calculate and return the micron-to-pixel ratio.
    """
    while True:
        try:
            micron_scale = float(input("Enter the known distance in microns (µm): "))
            pixels = float(input("Enter the corresponding distance in pixels: "))
            if pixels <= 0:
                print("Pixel distance must be positive. Try again.")
                continue
            ratio = micron_scale / pixels
            print(f"Micron-to-pixel ratio calculated: {ratio:.6f} µm/pixel\n")
            return ratio
        except ValueError:
            print("Invalid input. Please enter numeric values.")

def segment_sem_particles_tuned(
    image_path,
    microns_per_pixel,
    min_area_um2=0.50,
    max_area_um2=5.0
):
    """
    Tuned SEM nanoparticle segmentation to approach ImageJ results (~630 particles, 1.63 µm²).
    Returns an overlay image with detected nanoparticles outlined in red.
    """

    # Load image and crop bottom 10%
    img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise ValueError(f"Failed to load image: {image_path}")
    h, w = img.shape
    img = img[: int(h * 0.90), :]

    # Enhance contrast
    clahe = cv2.createCLAHE(clipLimit=2.5, tileGridSize=(8, 8))
    img_enh = clahe.apply(img)

    # Slight Gaussian blur (retain small particles, reduce noise)
    blurred = cv2.GaussianBlur(img_enh, (3, 3), 0)

    # Adaptive threshold
    binary = cv2.adaptiveThreshold(
        blurred, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
        cv2.THRESH_BINARY_INV, 15, 2
    )

    # Minimal morphology to separate touching particles without shrinking
    kernel = np.ones((3, 3), np.uint8)
    morph = cv2.morphologyEx(binary, cv2.MORPH_OPEN, kernel, iterations=1)
    morph = cv2.morphologyEx(morph, cv2.MORPH_CLOSE, kernel, iterations=1)

    # Distance transform + watershed for cluster splitting
    dist = cv2.distanceTransform(morph, cv2.DIST_L2, 5)
    _, sure_fg = cv2.threshold(dist, 0.05 * dist.max(), 255, 0)
    sure_fg = np.uint8(sure_fg)
    unknown = cv2.subtract(morph, sure_fg)

    # Marker labeling
    _, markers = cv2.connectedComponents(sure_fg)
    markers = markers + 1
    markers[unknown == 255] = 0

    # Watershed
    color_img = cv2.cvtColor(img_enh, cv2.COLOR_GRAY2BGR)
    cv2.watershed(color_img, markers)
    mask = np.uint8(markers > 1) * 255

    # Find contours
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Convert area thresholds
    min_area_px = min_area_um2 / (microns_per_pixel ** 2)
    max_area_px = max_area_um2 / (microns_per_pixel ** 2)

    valid_contours = []
    areas_um2 = []

    for c in contours:
        area_px = cv2.contourArea(c)
        if min_area_px <= area_px <= max_area_px:
            areas_um2.append(area_px * (microns_per_pixel ** 2))
            valid_contours.append(c)

    n_particles = len(areas_um2)
    avg_area = np.mean(areas_um2) if n_particles > 0 else 0.0

    print(f"Detected nanoparticles: {n_particles}")
    print(f"Average particle area: {avg_area:.4f} µm²")

    # Draw red contours on a copy of the grayscale image
    overlay = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    cv2.drawContours(overlay, valid_contours, -1, (0, 0, 255), 1)

    return overlay, areas_um2, avg_area


# Ask the user for micron-to-pixel ratio
microns_per_pixel = get_micron_to_pixel_ratio()

# Example usage:
image_input_path = "C:/Users/D00456326/Desktop/Research Photos/Nanoparticles/WCOB2.jpg"

overlay, areas, avg = segment_sem_particles_tuned(
    image_input_path,
    microns_per_pixel=microns_per_pixel,
    min_area_um2=0.50,
    max_area_um2=7.5
)

# Display the result
cv2.imshow("Detected Nanoparticles", overlay)
cv2.waitKey(0)
cv2.destroyAllWindows()
