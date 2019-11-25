# SEN1-ARD-GAMMA

This is a python module and terminal commands for processing Sentinel-1 GRD data products using the Gamma software. 

If you have Gamma installed locally this module can be run with the following command:

    sen1_grd_ard.py -i ./S1A_IW_GRDH_1SDV_20191124T063118_20191124T063143_030049_036E66_30F8.SAFE/ -o ./output/ -t ./tmp/ -d DTM_UK.kea -r 20 -f GTIFF --intimgs
    
If you have Gamma in a singularity image then you'll need to specify the following environmental variable

    export S1ARD_PAP_CMD="singularity exec --bind /data:/data /data/gamma_sfw.sif"
    sen1_grd_ard.py -i ./S1A_IW_GRDH_1SDV_20191124T063118_20191124T063143_030049_036E66_30F8.SAFE/ -o ./output/ -t ./tmp/ -d DTM_UK.kea -r 20 -f GTIFF --intimgs

With Docker you might also need to specify a local mount path which will be replaced when excuted:

    export S1ARD_PAP_CMD="docker run -i -t -v /data/sen1_analysis:/data gamma-sfw"
    export S1ARD_PAP_PATH="/data/sen1_analysis:/data"

    python sen1_grd_ard.py -i /data/sen1_analysis/inputs/S1A_IW_GRDH_1SDV_20181127T175729_20181127T175754_024777_02BA1D_432E.SAFE -o /data/sen1_analysis/outputs -t /data/sen1_analysis/tmp -d /data/sen1_analysis/srtm_global_mosaic_1arc_v3.tif -r 20  --intimgs

