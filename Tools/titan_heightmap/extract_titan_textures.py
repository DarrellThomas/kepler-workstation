#!/usr/bin/env python3
"""
Extract Titan surface textures from Cassini SAR and ISS mosaics.

Crops the global GeoTIFF mosaics to the Huygens landing site region,
generates a normal map from SAR backscatter, and exports 4096x4096 PNGs
ready for UE5 material import.

Input:  Cassini SAR global mosaic (351 m/px GeoTIFF)
        Cassini ISS global mosaic (4 km/px or 450 m/px GeoTIFF)
Output: T_TitanSAR_D.png  (SAR diffuse/detail texture)
        T_TitanISS_D.png  (ISS optical color texture)
        T_TitanSAR_N.png  (Normal map from SAR)

Usage:
    python3 extract_titan_textures.py [--size 4096] [--extent-km 719]
"""

import argparse
import json
import os
import sys
import numpy as np
from pathlib import Path

try:
    import rasterio
    from rasterio.windows import from_bounds
    from rasterio.transform import array_bounds
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False

try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

# Huygens landing site
HUYGENS_LAT = -10.573
HUYGENS_LON_W = 192.335
HUYGENS_LON_E = 360.0 - HUYGENS_LON_W  # 167.665°E

# Titan radius
TITAN_R_KM = 2575.0


def extract_geotiff_region(tiff_path, center_lon_e, center_lat, extent_km, output_size=4096):
    """Extract a region from a GeoTIFF and return as numpy array."""
    print(f"  Reading {os.path.basename(tiff_path)}...")

    with rasterio.open(tiff_path) as src:
        print(f"    CRS: {src.crs}")
        print(f"    Size: {src.width}x{src.height}")
        print(f"    Bounds: {src.bounds}")
        print(f"    Resolution: {src.res}")

        # Determine coordinate system
        # Most Titan mosaics use simple cylindrical (lon/lat in degrees)
        # or pixel coordinates. Check bounds to understand.
        bounds = src.bounds
        b_width = bounds.right - bounds.left
        b_height = bounds.top - bounds.bottom

        # If bounds are in degrees (typical for planetary GeoTIFFs)
        if b_width <= 360 and b_height <= 180:
            print(f"    Coordinate system: degrees (lon/lat)")
            # Convert extent from km to degrees at Titan's equator
            deg_per_km = 360.0 / (2 * np.pi * TITAN_R_KM)
            half_extent_deg = (extent_km / 2) * deg_per_km

            # Adjust for latitude (longitude degrees shrink with cos(lat))
            half_extent_lon = half_extent_deg / max(np.cos(np.radians(center_lat)), 0.1)
            half_extent_lat = half_extent_deg

            lon_min = center_lon_e - half_extent_lon
            lon_max = center_lon_e + half_extent_lon
            lat_min = center_lat - half_extent_lat
            lat_max = center_lat + half_extent_lat

            print(f"    Crop window: lon [{lon_min:.2f}, {lon_max:.2f}], lat [{lat_min:.2f}, {lat_max:.2f}]")

            # Clamp to data bounds
            lon_min = max(lon_min, bounds.left)
            lon_max = min(lon_max, bounds.right)
            lat_min = max(lat_min, bounds.bottom)
            lat_max = min(lat_max, bounds.top)

            try:
                window = from_bounds(lon_min, lat_min, lon_max, lat_max, src.transform)
                # Ensure integer window
                window = window.round_offsets().round_lengths()
                data = src.read(1, window=window)
                print(f"    Extracted: {data.shape[0]}x{data.shape[1]} pixels")
            except Exception as e:
                print(f"    Window extraction failed: {e}")
                print(f"    Falling back to center crop...")
                data = _center_crop(src, output_size)
        else:
            # Pixel coordinates or projected — just center crop
            print(f"    Coordinate system: projected/pixel (bounds: {b_width:.0f}x{b_height:.0f})")
            data = _center_crop(src, output_size)

    return data


def _center_crop(src, target_size):
    """Fallback: crop from center of the raster."""
    h, w = src.height, src.width
    # Read a square crop from near the Huygens region
    # Approximate: Huygens at ~47% from left, ~56% from top in global mosaic
    cx = int(w * (HUYGENS_LON_E / 360.0))
    cy = int(h * (90.0 - HUYGENS_LAT) / 180.0)

    half = min(target_size, min(h, w) // 2)
    y0 = max(0, cy - half)
    x0 = max(0, cx - half)
    y1 = min(h, y0 + 2 * half)
    x1 = min(w, x0 + 2 * half)

    window = rasterio.windows.Window(x0, y0, x1 - x0, y1 - y0)
    data = src.read(1, window=window)
    print(f"    Center crop: ({x0},{y0}) to ({x1},{y1}), {data.shape}")
    return data


def resize_array(data, target_size):
    """Resize 2D array to target_size x target_size using PIL."""
    if not HAS_PIL:
        # Numpy fallback: simple nearest-neighbor
        h, w = data.shape
        ys = np.linspace(0, h - 1, target_size).astype(int)
        xs = np.linspace(0, w - 1, target_size).astype(int)
        return data[np.ix_(ys, xs)]

    img = Image.fromarray(data)
    img = img.resize((target_size, target_size), Image.LANCZOS)
    return np.array(img)


def normalize_to_uint8(data):
    """Normalize data to 0-255 uint8."""
    data = data.astype(np.float64)
    # Handle NaN/nodata
    valid = np.isfinite(data) & (data > 0)  # Many mosaics use 0 as nodata
    if valid.sum() < 100:
        valid = np.isfinite(data)

    if valid.sum() == 0:
        return np.zeros(data.shape, dtype=np.uint8)

    vmin = np.percentile(data[valid], 2)
    vmax = np.percentile(data[valid], 98)

    if vmax - vmin < 1:
        vmin = data[valid].min()
        vmax = data[valid].max()

    clipped = np.clip(data, vmin, vmax)
    normalized = ((clipped - vmin) / (vmax - vmin) * 255).astype(np.uint8)
    return normalized


def generate_normal_map(data, strength=2.0):
    """Generate a normal map from elevation/intensity data using Sobel filter."""
    # Sobel kernels
    data_f = data.astype(np.float64)

    # Pad for edge handling
    padded = np.pad(data_f, 1, mode='edge')

    # Sobel X (horizontal gradient)
    gx = (padded[:-2, 2:] - padded[:-2, :-2] +
          2 * padded[1:-1, 2:] - 2 * padded[1:-1, :-2] +
          padded[2:, 2:] - padded[2:, :-2])

    # Sobel Y (vertical gradient)
    gy = (padded[2:, :-2] - padded[:-2, :-2] +
          2 * padded[2:, 1:-1] - 2 * padded[:-2, 1:-1] +
          padded[2:, 2:] - padded[:-2, 2:])

    gx *= strength
    gy *= strength

    # Compute normal vectors
    gz = np.ones_like(gx)
    length = np.sqrt(gx**2 + gy**2 + gz**2)
    length = np.where(length == 0, 1, length)

    nx = (-gx / length * 0.5 + 0.5) * 255
    ny = (-gy / length * 0.5 + 0.5) * 255
    nz = (gz / length * 0.5 + 0.5) * 255

    normal_map = np.stack([nx, ny, nz], axis=-1).astype(np.uint8)
    return normal_map


def apply_titan_colorize(sar_data):
    """Apply Titan's characteristic orange-brown color to SAR grayscale data.

    The SAR mosaic is grayscale (radar backscatter). We colorize it to match
    Huygens DISR imagery: dark tholin-orange plains with lighter ice features.
    """
    # Normalize to 0-1
    data_f = sar_data.astype(np.float64) / 255.0

    # Titan color palette (from Huygens POSTER_F)
    # Dark areas (low backscatter = smooth tholins): RGB ~(89, 46, 20) → (0.35, 0.18, 0.08)
    # Bright areas (high backscatter = rough ice): RGB ~(140, 115, 90) → (0.55, 0.45, 0.35)
    dark_rgb = np.array([89, 46, 20], dtype=np.float64)
    bright_rgb = np.array([145, 120, 90], dtype=np.float64)

    # Interpolate between dark and bright based on SAR intensity
    r = (dark_rgb[0] + data_f * (bright_rgb[0] - dark_rgb[0])).astype(np.uint8)
    g = (dark_rgb[1] + data_f * (bright_rgb[1] - dark_rgb[1])).astype(np.uint8)
    b = (dark_rgb[2] + data_f * (bright_rgb[2] - dark_rgb[2])).astype(np.uint8)

    return np.stack([r, g, b], axis=-1)


def main():
    parser = argparse.ArgumentParser(description='Extract Titan surface textures')
    parser.add_argument('--data-dir', default='../../PlanetaryModels/data/titan_photos',
                        help='Path to titan_photos directory')
    parser.add_argument('--output-dir', default='output',
                        help='Output directory')
    parser.add_argument('--size', type=int, default=4096,
                        help='Output texture resolution (default: 4096)')
    parser.add_argument('--extent-km', type=float, default=719.0,
                        help='Spatial extent in km (should match heightmap)')
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    data_dir = Path(args.data_dir)
    if not data_dir.is_absolute():
        data_dir = script_dir / data_dir

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = script_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if not HAS_RASTERIO:
        print("ERROR: rasterio required. pip install rasterio")
        sys.exit(1)

    print("=" * 60)
    print("  TITAN TEXTURE EXTRACTOR")
    print("=" * 60)
    print(f"\n  Center: Huygens ({HUYGENS_LAT}°, {HUYGENS_LON_E}°E)")
    print(f"  Extent: {args.extent_km:.0f} km")
    print(f"  Output: {args.size}×{args.size}")
    print()

    # ── SAR Mosaic (primary detail texture) ──
    sar_path = data_dir / 'cassini_radar' / 'Titan_SAR_Global_351m.tif'
    if sar_path.exists():
        print("━━━ SAR TEXTURE ━━━")
        sar_raw = extract_geotiff_region(str(sar_path), HUYGENS_LON_E, HUYGENS_LAT,
                                          args.extent_km, args.size)
        sar_uint8 = normalize_to_uint8(sar_raw)
        sar_resized = resize_array(sar_uint8, args.size)

        # Save grayscale SAR
        if HAS_PIL:
            Image.fromarray(sar_resized).save(str(output_dir / 'T_TitanSAR_Gray.png'))
            print(f"  Saved: T_TitanSAR_Gray.png ({args.size}×{args.size})")

        # Colorized version (Titan orange palette)
        sar_color = apply_titan_colorize(sar_resized)
        if HAS_PIL:
            Image.fromarray(sar_color).save(str(output_dir / 'T_TitanSAR_D.png'))
            print(f"  Saved: T_TitanSAR_D.png ({args.size}×{args.size}, colorized)")

        # Normal map from SAR
        print("  Generating normal map from SAR...")
        normal = generate_normal_map(sar_resized, strength=2.0)
        normal_resized = normal  # already at correct size
        if HAS_PIL:
            Image.fromarray(normal_resized).save(str(output_dir / 'T_TitanSAR_N.png'))
            print(f"  Saved: T_TitanSAR_N.png ({args.size}×{args.size})")
    else:
        print(f"  WARNING: SAR mosaic not found at {sar_path}")
        sar_resized = None

    print()

    # ── ISS Optical Mosaic (color reference) ──
    # Try high-res first, fall back to low-res
    iss_hires = data_dir / 'cassini_iss' / 'Titan_ISS_Globe_450m.tif'
    iss_lores = data_dir / 'cassini_iss' / 'Titan_ISS_Global_4km.tif'
    iss_path = iss_hires if iss_hires.exists() else iss_lores

    if iss_path.exists():
        print("━━━ ISS OPTICAL TEXTURE ━━━")
        iss_raw = extract_geotiff_region(str(iss_path), HUYGENS_LON_E, HUYGENS_LAT,
                                          args.extent_km, args.size)
        iss_uint8 = normalize_to_uint8(iss_raw)
        iss_resized = resize_array(iss_uint8, args.size)

        # ISS is grayscale (938nm methane window). Colorize with Titan palette.
        iss_color = apply_titan_colorize(iss_resized)
        if HAS_PIL:
            Image.fromarray(iss_color).save(str(output_dir / 'T_TitanISS_D.png'))
            print(f"  Saved: T_TitanISS_D.png ({args.size}×{args.size}, colorized)")
    else:
        print(f"  WARNING: ISS mosaic not found at {iss_path}")

    print()

    # ── Summary ──
    output_files = list(output_dir.glob('T_Titan*.png'))
    print("━━━ OUTPUT SUMMARY ━━━")
    for f in sorted(output_files):
        size_mb = f.stat().st_size / 1e6
        print(f"  {f.name}: {size_mb:.1f} MB")

    print(f"""
  ✓ DONE

  Textures ready for UE5 import:
    T_TitanSAR_D.png  — SAR diffuse (colorized Titan orange palette)
    T_TitanSAR_Gray.png — SAR grayscale (raw backscatter)
    T_TitanSAR_N.png  — Normal map (generated from SAR)
    T_TitanISS_D.png  — ISS optical (colorized)

  Color palette based on Huygens POSTER_F:
    Dark tholin plains: RGB(89, 46, 20)
    Bright ice features: RGB(145, 120, 90)
""")


if __name__ == '__main__':
    main()
