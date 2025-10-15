# IFC to Structured3D Converter

A tool for converting Industry Foundation Classes (IFC) building models to Structured3D format for visual localization research.

## Overview

This converter extracts geometric and semantic information from IFC files and transforms them into the Structured3D annotation format. It processes building elements (walls, doors, windows, spaces) and generates:

- **Junctions**: 3D coordinates of geometric vertices
- **Lines**: Edge connections between junctions
- **Planes**: Wall and opening surfaces with normals and offsets
- **Semantics**: Room and opening type annotations

The output format is compatible with the [SPVLoc](https://github.com/microsoft/spvloc) visual localization framework.

## Research

This work was presented at **Forum Bauinformatik 2025**. The full paper is available in the `paper/` directory.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/ifc-to-structured3d.git
cd ifc-to-structured3d

# Install dependencies
pip install -r requirements.txt
```

## Usage

1. Edit the file paths in `ifc_to_s3d.py`:
```python
ifc_file_path = "path/to/your/model.ifc"
output_json_path = "path/to/output/annotation_3d.json"
```

2. Run the converter:
```bash
python ifc_to_s3d.py
```

## Output Format

The converter generates a JSON file with the following structure:

```json
{
  "junctions": [...],      // 3D vertex coordinates
  "lines": [...],          // Edge definitions
  "planes": [...],         // Surface planes (walls, doors, windows)
  "semantics": [...]       // Room and opening annotations
}
```

## SPVLoc Modifications

The `spvloc_modifications/` directory contains custom modifications made to the original SPVLoc repository to support IFC-based data. These changes were necessary to adapt the visual localization pipeline for building information models.

See `spvloc_modifications/README.md` for details on what was modified and why.

## Dependencies

- `ifcopenshell` - IFC file parsing
- `numpy` - Numerical computations
- `scipy` - Spatial algorithms (ConvexHull)

## License

MIT License

## Citation

If you use this work, please cite our paper:

```bibtex
@inproceedings{soultana2025ifc,
  title={IFC to Structured3D Conversion for Visual Localization},
  author={Your Name},
  booktitle={Forum Bauinformatik 2025},
  year={2025}
}
```
