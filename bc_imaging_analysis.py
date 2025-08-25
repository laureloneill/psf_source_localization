"""BlackCAT Imaging Analysis, requires xraysky to be installed for imaging analysis to be run"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from pathlib import Path
import json

from xhcd.blackcat import calibrationAnalysis
from xraysky.scripts import events_imaging

def imaging_analysis(
    analysis_dir,
    mode='Full',
    th1_inint = 50,
    source = 'Ti',
):

    basedir = Path(analysis_dir)
    data_dir = Path(
        "/storage/group/adf15/default/Astro_X_Ray/HCD_Data/Speedster_data/Speedster550_data/"
    )
    
    
    det_sns = [23200, 23196, 23056, 23197]

    eventlist_dict = {}
    gainmaps_dict = {}
    badpix_masks = {}
    for sn in det_sns:

        eventlist_dict[sn] = np.sort(
            [f for f in basedir.glob(f"eventlist_s{sn}_th1_{th1_inint}*.fits*")]
        )

        if mode.lower() == 'event_driven':
            gainmaps_dict[sn] = (
                data_dir / f"s{sn}/gainmaps/Flight_gainmap_233K_25V_152Hz_Sparse_Mode.fits"
            )
        elif mode.lower() == 'event_driven_fast':
            gainmaps_dict[sn] = (
                data_dir / f"s{sn}/gainmaps/Flight_gainmap_233K_25V_700Hz_Sparse_Mode.fits"
            )
        else:
            gainmaps_dict[sn] = (
                data_dir / f"s{sn}/gainmaps/Flight_gainmap_233K_25V_Full_Frame.fits"
            )

        bp = fits.getdata(
            f"/storage/group/adf15/default/Astro_BlackCAT/BlackCAT_Calibration_Data/Badpixel_Analysis/Bad_Pixel_Masks/s{sn}/Flight_badpix_RTN_40e_233K_25V_Full_Frame.fits"
        )
        badpix_masks[sn] = np.asarray(bp,dtype=bool)

    outname = "imaging_analysis"

    th1_e = 200
    
    #print(eventlist_dict)
    
    analysis_args = {
        "eventlists": eventlist_dict,
        "gain_maps": gainmaps_dict,
        "badpix_maps": badpix_masks,
        "det_sn": det_sns,
        "det_pos_x": [1, 1, 0, 0],
        "det_pos_y": [0, 1, 0, 1],
        "run_dir": basedir,
        "th1": th1_e,
        "th2": 60,
        "source": source,
        "outname": outname
    }

    output_dir = basedir/f"Analysis/{outname}"
    
    analysis, analysis_res = calibrationAnalysis.runAnalysis(**analysis_args)

    imager = events_imaging.BC_Imaging()
    im_hdu = imager.evtlist2image(
        output_dir/f"combined_events_th1_{th1_e}.fits.gz",
        output_dir/f"image_reconstruction.fits.gz",
    )

    imager.plotSource(im_hdu, outname = "imaging", outdir = output_dir/"Figures")
    
    peaks = imager.imager.findpeaks(im_hdu.data.T)
    if len(peaks)==1:
        peakcoords = peaks[0]["radecpeak"]
        analysis_res['peakloc'] = list(peakcoords)
        analysis_res['peaksig'] = peaks[0]["signif_local"]

    with open(output_dir/'analysis_results.json','w') as f:
        json.dump(analysis_res,f)