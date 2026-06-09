#!/usr/bin/env python3
"""
Build UE5-ready heightmap from Cassini GTDR Titan topography data.

Input:  GTDR PDS3 IMG file (1440x1440, float32, equirectangular)
Output: .r16 heightmap (uint16 little-endian) for UE5 Landscape import
        metadata JSON with scale factors for UE5 configuration

Source: NASA PDS CO-SSA-RADAR-5-GTDR-V1.0 (Zebker et al., Stanford)

Usage:
    python3 build_titan_heightmap.py [--size 2017] [--region 128]
"""

import gzip
import struct
import json
import math
import argparse
import os
import sys
import numpy as np
from pathlib import Path

# GTDR file parameters (from PDS3 label)
GTDR_LINES = 1440
GTDR_SAMPLES = 1440
GTDR_SAMPLE_BYTES = 4  # float32
GTDR_RECORD_BYTES = 5760  # = 1440 * 4
GTDR_LABEL_RECORDS = 1
GTDR_HEADER_BYTES = GTDR_LABEL_RECORDS * GTDR_RECORD_BYTES  # 5760
GTDR_PIX_PER_DEG = 8.0
GTDR_KM_PER_PIX = 5.61777852985675
GTDR_MISSING = struct.unpack('<f', bytes.fromhex('FBFFFF7F'))[0]  # 16#FF7FFFFB# little-endian
GTDR_REF_RADIUS_KM = 2575.0

# N270 file covers: lon 180W-360W, lat 90N-90S
# POSITIVE_LONGITUDE_DIRECTION = WEST
# CENTER_LONGITUDE = 180W, SAMPLE_PROJECTION_OFFSET = 1439.5
# LINE_PROJECTION_OFFSET = 719.5

# Huygens landing site
HUYGENS_LAT = -10.573   # degrees (south)
HUYGENS_LON_W = 192.335  # degrees west


def read_gtdr(filepath):
    """Read a GTDR PDS3 IMG file (gzipped) and return elevation array in meters."""
    print(f"  Reading {filepath}...")

    open_fn = gzip.open if filepath.endswith('.gz') else open
    with open_fn(filepath, 'rb') as f:
        # Skip label record(s)
        f.read(GTDR_HEADER_BYTES)

        # Read image data
        raw = f.read(GTDR_LINES * GTDR_SAMPLES * GTDR_SAMPLE_BYTES)

    # Parse as float32 little-endian
    data = np.frombuffer(raw, dtype='<f4').reshape((GTDR_LINES, GTDR_SAMPLES))

    # Replace missing values with NaN
    missing_val = np.float32(GTDR_MISSING)
    # Also catch any very large/small sentinel values
    data = np.where(np.abs(data) > 1e6, np.nan, data)

    n_valid = np.count_nonzero(~np.isnan(data))
    n_total = data.size
    print(f"  Valid pixels: {n_valid}/{n_total} ({n_valid/n_total*100:.1f}%)")
    print(f"  Elevation range: {np.nanmin(data):.1f} to {np.nanmax(data):.1f} m")

    return data


def lonlat_to_pixel(lon_w, lat, offset_sample=1439.5, offset_line=719.5):
    """Convert west-longitude/latitude to pixel coordinates in the N270 tile.

    N270 tile: westernmost=360W, easternmost=180W
    Sample increases with increasing west longitude from 180W.
    """
    sample = (lon_w - 180.0) * GTDR_PIX_PER_DEG
    line = (90.0 - lat) * GTDR_PIX_PER_DEG
    return int(round(line)), int(round(sample))


def extract_region(data, center_lat, center_lon_w, half_size_pixels):
    """Extract a square region centered on a coordinate."""
    cy, cx = lonlat_to_pixel(center_lon_w, center_lat)
    print(f"  Huygens pixel coords: line={cy}, sample={cx}")

    y0 = max(0, cy - half_size_pixels)
    y1 = min(GTDR_LINES, cy + half_size_pixels)
    x0 = max(0, cx - half_size_pixels)
    x1 = min(GTDR_SAMPLES, cx + half_size_pixels)

    region = data[y0:y1, x0:x1].copy()

    # Fill NaN gaps with interpolation
    nan_mask = np.isnan(region)
    n_nan = np.count_nonzero(nan_mask)
    if n_nan > 0:
        print(f"  Filling {n_nan} NaN pixels with mean elevation...")
        mean_val = np.nanmean(region)
        region[nan_mask] = mean_val

    extent_km = region.shape[0] * GTDR_KM_PER_PIX
    print(f"  Region: {region.shape[0]}x{region.shape[1]} pixels = {extent_km:.0f}x{extent_km:.0f} km")
    print(f"  Elevation range: {region.min():.1f} to {region.max():.1f} m")

    return region, (y0, x0, y1, x1)


def upsample_bicubic(data, target_size):
    """Bicubic upsample using scipy if available, else bilinear with numpy."""
    try:
        from scipy.ndimage import zoom
        factor = target_size / data.shape[0]
        result = zoom(data, factor, order=3)  # bicubic
        # Ensure exact target size
        result = result[:target_size, :target_size]
        print(f"  Upsampled {data.shape[0]}x{data.shape[1]} → {result.shape[0]}x{result.shape[1]} (bicubic)")
        return result
    except ImportError:
        # Fallback: numpy bilinear
        print("  scipy not available, using bilinear interpolation...")
        h, w = data.shape
        new_h, new_w = target_size, target_size
        y_ratio = h / new_h
        x_ratio = w / new_w

        result = np.zeros((new_h, new_w), dtype=np.float64)
        for i in range(new_h):
            for j in range(new_w):
                y = i * y_ratio
                x = j * x_ratio
                y0, x0 = int(y), int(x)
                y1, x1 = min(y0 + 1, h - 1), min(x0 + 1, w - 1)
                fy, fx = y - y0, x - x0
                result[i, j] = (data[y0, x0] * (1 - fy) * (1 - fx) +
                                data[y1, x0] * fy * (1 - fx) +
                                data[y0, x1] * (1 - fy) * fx +
                                data[y1, x1] * fy * fx)
        print(f"  Upsampled {h}x{w} → {new_h}x{new_w} (bilinear)")
        return result


def add_procedural_detail(data, amplitude_m=30.0, octaves=6, seed=42):
    """Add fractal noise detail below GTDR resolution.

    Uses multi-octave Perlin-like noise to add plausible terrain detail
    at scales smaller than the 5.62 km GTDR pixel size.
    """
    np.random.seed(seed)
    h, w = data.shape
    noise = np.zeros((h, w))

    for octave in range(octaves):
        freq = 2 ** (octave + 2)
        amp = amplitude_m / (2 ** octave)

        # Generate smooth noise at this frequency
        small = np.random.randn(freq, freq)
        # Upsample to full resolution
        from scipy.ndimage import zoom as sz
        try:
            layer = sz(small, h / freq, order=3)[:h, :w]
        except Exception:
            # Fallback: tile and smooth
            layer = np.tile(small, (h // freq + 1, w // freq + 1))[:h, :w]

        noise += layer * amp

    print(f"  Added procedural detail: {octaves} octaves, max amplitude ±{amplitude_m:.0f} m")
    return data + noise


def export_r16(data, filepath, metadata_path=None):
    """Export as .r16 (uint16 little-endian) for UE5 Landscape import."""
    elev_min = data.min()
    elev_max = data.max()
    elev_range = elev_max - elev_min

    # Map elevation to uint16 (0-65535)
    normalized = (data - elev_min) / elev_range
    uint16_data = (normalized * 65535).astype(np.uint16)

    # Write raw uint16 little-endian
    uint16_data.tofile(filepath)

    size_mb = os.path.getsize(filepath) / 1e6
    print(f"  Exported {filepath}: {data.shape[0]}x{data.shape[1]}, {size_mb:.1f} MB")

    # Metadata for UE5 import
    metadata = {
        'width': int(data.shape[1]),
        'height': int(data.shape[0]),
        'format': 'uint16_le',
        'elevation_min_m': float(elev_min),
        'elevation_max_m': float(elev_max),
        'elevation_range_m': float(elev_range),
        'ue5_z_scale': float(elev_range * 100 / 65535),  # cm per uint16 unit
        'center_lat': HUYGENS_LAT,
        'center_lon_w': HUYGENS_LON_W,
        'gtdr_km_per_pixel': GTDR_KM_PER_PIX,
        'output_m_per_pixel': float(elev_range),  # will be overridden below
        'titan_ref_radius_km': GTDR_REF_RADIUS_KM,
    }

    if metadata_path:
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        print(f"  Metadata: {metadata_path}")

    return metadata


def main():
    parser = argparse.ArgumentParser(description='Build Titan heightmap from GTDR data')
    parser.add_argument('--gtdr-dir', default='../../PlanetaryModels/data/titan_topography/gtdr-data',
                        help='Path to GTDR data directory')
    parser.add_argument('--output-dir', default='output',
                        help='Output directory for .r16 and metadata')
    parser.add_argument('--region', type=int, default=128,
                        help='Half-size of extraction region in GTDR pixels (default: 128 = 1436 km)')
    parser.add_argument('--size', type=int, default=2017,
                        help='Output heightmap resolution (default: 2017, valid UE5 size)')
    parser.add_argument('--detail-amplitude', type=float, default=25.0,
                        help='Procedural detail max amplitude in meters (default: 25)')
    parser.add_argument('--no-detail', action='store_true',
                        help='Skip procedural detail (output only upsampled GTDR)')
    args = parser.parse_args()

    # Resolve paths
    script_dir = Path(__file__).parent
    gtdr_dir = Path(args.gtdr_dir)
    if not gtdr_dir.is_absolute():
        gtdr_dir = script_dir / gtdr_dir

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = script_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find the interpolated debiased N270 file (contains Huygens site)
    gtdr_file = gtdr_dir / 'GTIED00N270_T077_V01.IMG.gz'
    if not gtdr_file.exists():
        # Try uncompressed
        gtdr_file = gtdr_dir / 'GTIED00N270_T077_V01.IMG'
    if not gtdr_file.exists():
        print(f"ERROR: GTDR file not found: {gtdr_file}")
        print(f"  Expected in: {gtdr_dir}")
        sys.exit(1)

    print("=" * 60)
    print("  TITAN HEIGHTMAP BUILDER")
    print("=" * 60)
    print(f"\n  GTDR file: {gtdr_file}")
    print(f"  Huygens site: {HUYGENS_LAT}°, {HUYGENS_LON_W}°W")
    print(f"  Region: {args.region * 2} × {args.region * 2} GTDR pixels")
    print(f"  Output size: {args.size} × {args.size}")
    print()

    # Step 1: Read GTDR
    data = read_gtdr(str(gtdr_file))

    # Step 2: Extract region around Huygens
    print(f"\n  Extracting region around Huygens...")
    region, bounds = extract_region(data, HUYGENS_LAT, HUYGENS_LON_W, args.region)

    # Compute spatial extent
    extent_km = region.shape[0] * GTDR_KM_PER_PIX
    output_m_per_pixel = (extent_km * 1000) / args.size

    # Step 3: Upsample
    print(f"\n  Upsampling to {args.size}x{args.size}...")
    upsampled = upsample_bicubic(region, args.size)

    # Step 4: Add procedural detail
    if not args.no_detail:
        print(f"\n  Adding procedural terrain detail...")
        try:
            upsampled = add_procedural_detail(upsampled, amplitude_m=args.detail_amplitude)
        except ImportError:
            print("  WARNING: scipy needed for procedural detail, skipping")

    # Step 5: Export
    print(f"\n  Exporting...")
    r16_path = output_dir / 'titan_huygens_heightmap.r16'
    meta_path = output_dir / 'titan_huygens_heightmap.json'
    metadata = export_r16(upsampled, str(r16_path), str(meta_path))

    # Update metadata with spatial info
    metadata['output_m_per_pixel'] = output_m_per_pixel
    metadata['extent_km'] = extent_km
    metadata['ue5_xy_scale_cm_per_pixel'] = output_m_per_pixel * 100
    with open(str(meta_path), 'w') as f:
        json.dump(metadata, f, indent=2)

    # Also export a preview PNG
    try:
        from PIL import Image
        preview = ((upsampled - upsampled.min()) / (upsampled.max() - upsampled.min()) * 255).astype(np.uint8)
        Image.fromarray(preview).save(str(output_dir / 'titan_huygens_preview.png'))
        print(f"  Preview: {output_dir / 'titan_huygens_preview.png'}")
    except ImportError:
        print("  (PIL not available, skipping PNG preview)")

    print(f"""
  ✓ DONE

  Heightmap: {r16_path}
  Metadata:  {meta_path}

  UE5 Import Settings:
    File:           {r16_path.name}
    Resolution:     {args.size} × {args.size}
    Scale X/Y:      {output_m_per_pixel * 100:.1f} cm/pixel ({output_m_per_pixel:.1f} m/pixel)
    Scale Z:        {metadata['ue5_z_scale']:.2f} cm/unit (elevation range: {metadata['elevation_range_m']:.0f} m)
    Extent:         {extent_km:.0f} × {extent_km:.0f} km
    Center:         Huygens landing site ({HUYGENS_LAT}°, {HUYGENS_LON_W}°W)
""")


if __name__ == '__main__':
    main()
