# spe_parser
A simple library to load SPE file for spectroscopy


# How to install

If `git` is installed in your environment, 
```bash
pip install git+https://github.com/fujiisoup/spe_parser.git
```
to install `spe_parser`.

If not, download the package, and do the following
```bash
python setup.py install
```


# Basic usage
```python
import spe_parser
spe_parser.np_open('path/to/spe_file.spe')
```

If `xarray` is installed in your device, you can do
```python
spe_parser.xr_open('path/to/spe_file.spe')
```
which gives more complete information in the file.

The metadata can be obtained by
```python
data = spe_parser.xr_open('path/to/spe_file.spe', attributes='all')
```
then, the metadata will be stored in the `xr.DataArray`. You can see the contents by
```python
for key, item in data.attrs.items():
    print(key, item)
```

# Note
This package is still under the development. Any API may be the subject of the future change.
