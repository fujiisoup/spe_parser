"""
This module imports a Princeton Instruments LightField (SPE 3.0) file into a python environment.

Some part of the source code is bollowed from 
spe2py: https://github.com/ashirsch/spe2py.git
"""
from datetime import datetime
import numpy as np
import untangle
from io import StringIO


def robust_float(s):
    try:
        if '.' in s:
            return float(s)
        else:
            return int(s)
    except ValueError:
        return s


def np_open(spe_file):
    ''' Convert spe_file to np.ndarray
    All the data will be contained as the Coordinates or attrs

    Parameters
    ----------
    spe_file: filename
    attributes: string or a list of strings
        'default': attaches default attributes
        'all': attaches all the attributes
        list of strings: attaches the corresponding attributes
        dict: a mapping from original key to new key

    Returns
    -------
    array: np.ndarray
    info: dict of metadata
    '''
    data = SpeFile(spe_file)
    # TODO store metadata into a dict
    return data.data, []
    

def xr_open(spe_file, attributes='default'):
    """ Convert spe_file to xarray.DataArray
    All the data will be contained as the Coordinates or attrs

    Parameters
    ----------
    spe_file: filename
    attributes: string or a list of strings
        'default': attaches default attributes
        'all': attaches all the attributes
        list of strings: attaches the corresponding attributes
        dict: a mapping from original key to new key

    Returns
    -------
    xr.DataArray
    """
    try:
        import xarray as xr
    except ImportError:
        raise ImportError("To use `xr_open`, `xarray` needs to be installed to the system.")

    spefile = SpeFile(spe_file)
    coords = {}
    if spefile.nroi == len(spefile.data[0]):
        coords['x_original'] = [
            np.arange(0, int(spefile.roi[i]['width']), int(spefile.roi[i]['xBinning']))
            + int(spefile.roi[i]['x'])
            for i in range(spefile.nroi)
        ]
        xsizes = [len(x) for x in coords['x_original']]

        coords['y_original'] = [
            np.arange(0, int(spefile.roi[i]['height']), int(spefile.roi[i]['yBinning']))
            + int(spefile.roi[i]['y'])
            for i in range(spefile.nroi)
        ]
        ysizes = [len(y) for y in coords['y_original']]

    if spefile.wavelength is not None and len(spefile.wavelength) == len(coords['x_original'][0]):
        coords['wavelength'] = 'x', spefile.wavelength

    all_attrs = spefile.extract_all_metadata()
    # date
    date = all_attrs['SpeFormat.DataHistories.DataHistory.created']
    coords['date'] = np.datetime64(date[:-6])
    coords['date_timezone'] = {'-': -1, '+': 1}[date[-6]] * (
        np.timedelta64(date[-5:-3], 'h') + np.timedelta64(date[-2:], 'm')
    )

    attrs = {}
    if attributes == 'all':
        attrs = all_attrs
    elif attributes == 'default':
        key_pairs = {
            'SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera.Gating.RepetitiveGate.delay': 'gate_delay',
            'SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera.Gating.RepetitiveGate.width': 'gate_width'
        }
        for key, key_name in key_pairs.items():
            if key in all_attrs:
                attrs[key] = all_attrs[key]
    elif isinstance(attributes, (list, tuple)):
        for key in attributes:
            attrs[key] = all_attrs[key]
    elif isinstance(attributes, dict):
        for key, newkey in attributes.items():
            if key in all_attrs:
                attrs[newkey] = all_attrs[key]

    attrs['light_field_version'] = spefile.header_version

    if len(np.unique(xsizes)) == 1 and len(np.unique(ysizes)) == 1:
        coords['roi'] = np.arange(spefile.nroi)
        coords['x_original'] = ('roi', 'x'), coords['x_original']
        coords['y_original'] = ('roi', 'y'), coords['y_original']
        data = np.array(spefile.data)
        
        if data.shape[-1] != np.array(coords['x_original'][-1]).shape[-1]:
            coords['x_original'] = ('roi', 'x'), [[0]] * data.shape[-1]
        if data.shape[-2] != np.array(coords['y_original'][-1]).shape[-1]:
            coords['y_original'] = ('roi', 'y'), np.zeros((1, data.shape[-2]))
        return xr.DataArray(
            data, dims=['time', 'roi', 'y', 'x'],
            attrs=attrs, coords=coords
        )

    da = []
    for i in range(spefile.nroi):
        coords1 = coords.copy()
        coords1['x_original'] = ('x', ), coords['x_original'][i]
        coords1['y_original'] = ('y', ), coords['y_original'][i]
        coords1['roi'] = i
        da.append(xr.DataArray(
            [spefile.data[t][i] for t in range(len(spefile.data))],
            dims=['time', 'y', 'x'], coords=coords1,
            attrs=attrs
        ))

    if len(np.unique(xsizes)) == 1:
        return xr.concat([d for d in da], dim='y')
    elif len(np.unique(ysizes)) == 1:
        return xr.concat([d for d in da], dim='x')
    else:
        return xr.concat([d.stack(xy=['x', 'y']).reset_index('xy') for d in da], dim='xy')


class SpeFile:
    """
    Most of this class is bollowed from 
    spe2py: https://github.com/ashirsch/spe2py.git
    """
    def __init__(self, filepath):
        assert isinstance(filepath, str), 'Filepath must be a single string'
        self.filepath = filepath
 
        with open(self.filepath, encoding="utf8") as file:
            self.header_version = read_at(file, 1992, 3, np.float32)[0]
            assert self.header_version >= 3.0, \
                'This version of spe2py cannot load filetype SPE v. %.1f' % self.header_version

            self.nframes = read_at(file, 1446, 2, np.uint16)[0]

            self.footer = self._read_footer(file)
            self.dtype = self._get_dtype(file)

            # Note: these methods depend on self.footer
            self.xdim, self.ydim = self._get_dims()
            self.roi, self.nroi = self._get_roi_info()
            self.wavelength = self._get_wavelength()

            self.xcoord, self.ycoord = self._get_coords()

            self.data, self.metadata, self.metanames = self._read_data(file)
        file.close()

    @staticmethod
    def _read_footer(file):
        """
        Loads and parses the source file's xml footer metadata to an 'untangle' object.
        """
        footer_pos = read_at(file, 678, 8, np.uint64)[0]

        file.seek(footer_pos)
        xmltext = file.read()

        parser = untangle.make_parser()
        sax_handler = untangle.Handler()
        parser.setContentHandler(sax_handler)

        parser.parse(StringIO(xmltext))

        loaded_footer = sax_handler.root

        return loaded_footer

    @staticmethod
    def _get_dtype(file):
        """
        Returns the numpy data type used to encode the image data by reading the numerical code in the binary header.
        Reference: Princeton Instruments File Specification pdf
        """
        dtype_code = read_at(file, 108, 2, np.uint16)[0]

        if dtype_code == 0:
            dtype = np.float32
        elif dtype_code == 1:
            dtype = np.int32
        elif dtype_code == 2:
            dtype = np.int16
        elif dtype_code == 3:
            dtype = np.uint16
        elif dtype_code == 8:
            dtype = np.uint32
        else:
            raise ValueError("Unrecognized data type code: %.2f. Value should be one of {0, 1, 2, 3, 8}" % dtype_code)

        return dtype

    def _get_meta_dtype(self):
        meta_types = []
        meta_names = []
        prev_item = None
        for item in dir(self.footer.SpeFormat.MetaFormat.MetaBlock):
            if item == 'TimeStamp' and prev_item != 'TimeStamp':  # Specify ExposureStarted vs. ExposureEnded
                for element in self.footer.SpeFormat.MetaFormat.MetaBlock.TimeStamp:
                    meta_names.append(element['event'])
                    meta_types.append(element['type'])
                prev_item = 'TimeStamp'
            elif item == 'GateTracking' and prev_item != 'GateTracking':  # Specify Delay vs. Width
                for element in self.footer.SpeFormat.MetaFormat.MetaBlock.GateTracking:
                    meta_names.append(element['component'])
                    meta_types.append(element['type'])
                prev_item = 'GateTracking'
            elif prev_item != item:  # All other metablock names only have one possible value
                meta_names.append(item)
                meta_types.append(getattr(self.footer.SpeFormat.MetaFormat.MetaBlock, item)['type'])
                prev_item = item

        for index, type_str in enumerate(meta_types):
            if type_str == 'Int64':
                meta_types[index] = np.int64
            else:
                meta_types[index] = np.float64

        return meta_types, meta_names

    def _get_roi_info(self):
        """
        Returns region of interest attributes and numbers of regions of interest
        """
        try:
            camerasettings = self.footer.SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera
            regionofinterest_selection = camerasettings.ReadoutControl.RegionsOfInterest.Selection
            regionofinterest = camerasettings.ReadoutControl.RegionsOfInterest.CustomRegions.RegionOfInterest
            
            if not isinstance(regionofinterest, list):
                regionofinterest = [regionofinterest]
                nroi = 1
            else:
                nroi = len(regionofinterest)

            if regionofinterest_selection.cdata == 'BinnedSensor':
                roi = []
                for regionofinterest1 in regionofinterest:
                    roi.append({})
                    for key in ['x', 'xBinning', 'width', 'y', 'yBinning', 'height']:
                        roi[-1][key] = regionofinterest1[key]
                    binned_sensor = camerasettings.ReadoutControl.RegionsOfInterest.BinnedSensor
                    roi[-1]['xBinning'] = int(binned_sensor.XBinning.cdata)
                    roi[-1]['yBinning'] = int(binned_sensor.YBinning.cdata)
            else:
                roi = regionofinterest
    
        except (AttributeError):
            print("XML Footer was not loaded prior to calling _get_roi_info")
            raise

        return roi, nroi

    def _get_wavelength(self):
        """
        Returns wavelength-to-pixel map as stored in XML footer
        """
        try:
            wavelength_string = StringIO(self.footer.SpeFormat.Calibrations.WavelengthMapping.Wavelength.cdata)
        except AttributeError:
            #print("XML Footer was not loaded prior to calling _get_wavelength")
            return
        except IndexError:
            #print("XML Footer does not contain Wavelength Mapping information")
            return

        wavelength = np.loadtxt(wavelength_string, delimiter=',')

        return wavelength

    def _get_dims(self):
        """
        Returns the x and y dimensions for each region as stored in the XML footer
        """
        xdim = [int(block["width"]) for block in self.footer.SpeFormat.DataFormat.DataBlock.DataBlock]
        ydim = [int(block["height"]) for block in self.footer.SpeFormat.DataFormat.DataBlock.DataBlock]

        return xdim, ydim

    def _get_coords(self):
        """
        Returns x and y pixel coordinates. Used in cases where xdim and ydim do not reflect image dimensions
        (e.g. files containing frames with multiple regions of interest)
        """
        xcoord = [[] for _ in range(0, self.nroi)]
        ycoord = [[] for _ in range(0, self.nroi)]

        for roi_ind in range(0, self.nroi):
            working_roi = self.roi[roi_ind]
            ystart = int(working_roi['y'])
            ybinning = int(working_roi['yBinning'])
            yheight = int(working_roi['height'])
            ycoord[roi_ind] = range(ystart, (ystart + yheight), ybinning)

        for roi_ind in range(0, self.nroi):
            working_roi = self.roi[roi_ind]
            xstart = int(working_roi['x'])
            xbinning = int(working_roi['xBinning'])
            xwidth = int(working_roi['width'])
            xcoord[roi_ind] = range(xstart, (xstart + xwidth), xbinning)

        return xcoord, ycoord

    def _read_data(self, file):
        """
        Loads raw image data into an nframes X nroi list of arrays.
        """
        file.seek(4100)

        frame_stride = int(self.footer.SpeFormat.DataFormat.DataBlock['stride'])
        frame_size = int(self.footer.SpeFormat.DataFormat.DataBlock['size'])
        metadata_size = frame_stride - frame_size
        if metadata_size != 0:
            metadata_dtypes, metadata_names = self._get_meta_dtype()
            metadata = np.zeros((self.nframes, len(metadata_dtypes)))
        else:
            metadata_dtypes, metadata_names = None, None
            metadata = None

        data = [[0 for _ in range(self.nroi)] for _ in range(self.nframes)]
        for frame in range(0, self.nframes):
            for region in range(0, self.nroi):
                if self.nroi > 1:
                    data_xdim = len(self.xcoord[region])
                    data_ydim = len(self.ycoord[region])
                else:
                    data_xdim = np.asarray(self.xdim[region], np.uint32)
                    data_ydim = np.asarray(self.ydim[region], np.uint32)
                data[frame][region] = np.fromfile(file, self.dtype, data_xdim * data_ydim).reshape(data_ydim, data_xdim)
            if metadata_dtypes is not None:
                for meta_block in range(len(metadata_dtypes)):
                    metadata[frame, meta_block] = np.fromfile(file, dtype=metadata_dtypes[meta_block], count=1)

        return data, metadata, metadata_names

    def xmltree(self, footer, ind=-1):
        """
        Prints the untangle footer object in tree form to easily view metadata fields. Ignores object elements that
        contain lists (e.g. ..Spectrometer.Turrets.Turret).
        """
        if dir(footer):
            ind += 1
            for item in dir(footer):
                if isinstance(getattr(footer, item), list):
                    continue
                else:
                    print(ind * ' -->', item)
                    self.xmltree(getattr(footer, item), ind)

    def extract_all_metadata(self):
        def extract(obj):
            dicts = {}
            keys = [k for k in obj.__dir__() if k[0] != '_']
            for k in keys:
                newobj = getattr(obj, k)
                if isinstance(newobj, list):
                    for i, newobj1 in enumerate(newobj):
                        for newkey, newitem in newobj1._attributes.items():
                            if newkey != 'id':
                                dicts[str(i) + '.' + newkey] = newitem
                elif isinstance(newobj, untangle.Element):
                    for newkey, newitem in newobj._attributes.items():
                        if newkey != 'id':
                            dicts[newkey] = newitem

                    newdict = extract(newobj)
                    for newkey, newitem in newdict.items():
                        dicts[k + '.' + newkey] = newitem
            return dicts
        
        return extract(self.footer)



def read_at(file, pos, size, ntype):
    """
    Reads SPE source file at specific byte position.
    Adapted from https://scipy.github.io/old-wiki/pages/Cookbook/Reading_SPE_files.html
    """
    file.seek(pos)
    return np.fromfile(file, ntype, size)

