from multiprocessing.sharedctypes import Value
import os
import sys
THIS_DIR = os.path.dirname(__file__)
sys.path.append(THIS_DIR + '/../spe_parser/')

# data directories that will be tested
PUBLICDATA_DIR = THIS_DIR + '/public_data/'
PRIVATEDATA_DIR = THIS_DIR + '/private_data/'

import numpy as np
import pytest
import spe_parser


filenames = [PUBLICDATA_DIR + f for f in os.listdir(PUBLICDATA_DIR) if f[-4:] in ['.spe', '.SPE']]
filenames += [PRIVATEDATA_DIR + f for f in os.listdir(PRIVATEDATA_DIR) if f[-4:] in ['.spe', '.SPE']]
print(filenames)

@pytest.mark.parametrize('filename', filenames)
def test_open(filename):
    data, info = spe_parser.np_open(filename)
    assert np.sum(np.isnan(data)) == 0

    # open reference data if exists
    reffile = filename[:-4] + '.npz'
    if os.path.exists(reffile):
        ref = np.load(reffile)
        ref = ref[list(ref.keys())[0]]
        assert np.allclose(ref, data)

    data = spe_parser.xr_open(filename)
    # testing attributes
    # no error
    data = spe_parser.xr_open(filename, attributes='all')
    
    # should have new key
    attrs = {
        'SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera.HardwareIO.AuxOutput.Gate.width':
        'gate_width'
    }
    data = spe_parser.xr_open(filename, attributes=attrs)
    assert 'gate_width' in data.attrs 
    
