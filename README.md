# BB86: HYCOM Model Utilities

## References:

- (BB86 For Dummies)[https://www.hycom.org/attachments/349_BB86_for_dummies.pdf]
- (Hycom for Dummies)[https://www.hycom.org/attachments/349_HYCOM_for_Dummies_2011.pdf]
- (HYCOM Documentation)[https://www.hycom.org/hycom/documentation]

## Usage

### Converting HYCOM archive ([a|b]binary) files to netcdf.

#### Generate Grid File

```bash
python3 /usr/local/BB86_PACKAGE/PYTHON/write_grid_bb86.py
```

#### Generate Depth File

```bash
python3 /usr/local/BB86_PACKAGE/PYTHON/write_depth_bb86.py
```

#### Convert Archives to Netcdf

```bash
csh /usr/local/hycom_example_scripts/archive/ncdf3z.com

```
