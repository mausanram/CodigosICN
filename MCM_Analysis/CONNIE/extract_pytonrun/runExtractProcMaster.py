import os
import glob
import sys
import subprocess
import re
import numpy as np
import math
import shlex
import argparse
import shutil

# Create the parser
my_parser = argparse.ArgumentParser(description='This program will run the exctraction CONNIE routine', epilog='Enjoy your analysis!')
my_parser.add_argument('-r','--runIDs', nargs='+',type=int, action='store', help='List of runID to analyse and produce .fits files. * Just write the runIDs number separeted by space.')
my_parser.add_argument('-v','--version', action='store', help='Version of extraction code.')
my_parser.add_argument('--exe', action='store', help='Path to extract.exe (default: ./extract.exe)')
my_parser.add_argument('-x','--channels', action='store', help='Comma-separated channels to pass to extract.exe (e.g. 4,7,9)')
my_parser.add_argument('--data-dir', required=True, action='store', help='Base data dir for images (e.g. /share/.../calib)')
my_parser.add_argument('--mask-dir', required=True, action='store', help='Base mask dir (e.g. /share/.../masksBL)')
my_parser.add_argument('--catalog-dir', action='store', help='Base output catalog dir (overrides default, default: ./catalogs)')
args = my_parser.parse_args()

for run in args.runIDs:
    
    # Data and mask paths must be provided (no old CONNIE defaults)
    DATAPATH = os.path.join(args.data_dir, str(run))
    MASKPATH = os.path.join(args.mask_dir, str(run))

    # Choose extract executable (default: ./extract.exe)
    if args.exe:
        EXE = args.exe
    else:
        EXE = os.path.join(os.getcwd(),'extract.exe')

    # Choose catalog output path (default: ./catalogs/<version>/<run> or ./catalogs/<run>)
    if args.catalog_dir:
        CATALOGPATH = os.path.join(args.catalog_dir, args.version if args.version else '', str(run))
    else:
        if args.version:
            CATALOGPATH=os.path.join(os.getcwd(),'catalogs',args.version,str(run))
        else:
            CATALOGPATH=os.path.join(os.getcwd(),'catalogs',str(run))
    seed='10'

    print('\033[92m'+'START PROCESSING RUNID %s ' % run + "\n \033[0m")

    # Creating individual catalogs by iterating masks
    try:
        os.makedirs(CATALOGPATH, exist_ok = True)
    except OSError as error:
        print("Directory '%s' can not be created" % CATALOGPATH)
        sys.exit()

    mask_files = sorted(glob.glob(os.path.join(MASKPATH,'*.fits')))
    if not mask_files:
        print('\033[93m'+'No .fits masks found in %s' % MASKPATH + '\033[0m')
        continue

    for path_mask in mask_files:
        mask_basename = os.path.basename(path_mask)

        # The cal filename always contains 'calib_equalized...' — extract that part from the mask name
        idx = mask_basename.find('calib_equalized')
        if idx != -1:
            cal_basename = mask_basename[idx:]
        else:
            # fallback: remove common prefix like 'BLmask_' or 'GBmask_'
            cal_basename = re.sub(r'^(?:BLmask_|GBmask_|mask_)','', mask_basename)

        cal = os.path.join(DATAPATH, cal_basename)
        base = os.path.splitext(mask_basename)[0]
        if args.version:
            catalog = os.path.join(CATALOGPATH, 'catalog_M2_{}_{}_{}s.root'.format(args.version, base, seed))
        else:
            catalog = os.path.join(CATALOGPATH, 'catalog_M2_{}_{}s.root'.format(base, seed))

        cmd_args = [EXE, cal, '-m', path_mask, '-r', str(run), '-o', catalog]
        if args.channels:
            cmd_args += ['-x', args.channels]

        print('\nRunning: {}'.format(' '.join(cmd_args)))
        proc=subprocess.Popen(cmd_args)
        proc.wait()

    print('Done! The catalogs will be merged into one: \n')

    # Merging catalogs
    merged_out = os.path.join(CATALOGPATH, ('catalog_M2_'+(args.version if args.version else 'run')+'_runID'+str(run)+'_'+seed+'s.root'))
    hadd_cmd = ['hadd','-f', merged_out] + sorted(glob.glob(os.path.join(CATALOGPATH,'catalog_M2_*_'+seed+'s.root')))
    if len(hadd_cmd) > 2:
        os.system(' '.join(hadd_cmd))
        print('Merged!')
    else:
        print('\033[93m'+'Not enough catalogs to merge in %s' % CATALOGPATH + '\033[0m')
    	
#     # Selection Energy Cut in the catalog (uses selEvents.exe in CONNIE tools)
#     SELE_PATH='/share/storage2/connie/ConnieCode/lowECatalog/v1.0'
#     merged_in = merged_out
#     low_out = os.path.join(CATALOGPATH, 'low_catalog_M2_'+(args.version if args.version else 'run')+'_runID'+str(run)+'_'+seed+'s.root')
#     sel_cmd = SELE_PATH + '/selEvents.exe -o ' + low_out + ' -s "E0<5000" -i ' + merged_in
#     os.system(sel_cmd)
    
#     # Perform Sigma x and y fit in the events
#    print('\n Performing fits to the events: \n')
#    pyPATH='/share/storage2/connie/ConnieCode/CONNIEtools'
#    command_line = 'python '+pyPATH+'/Catalog.py addbranches '+CATALOGPATH+'/low_catalog_M2_'+args.version+'_runID'+str(run)+'.root fit -o '+CATALOGPATH+'/sigma_low_catalog_M2_'+args.version+'_runID'+str(run)+'.root'
#    os.system(command_line)
    
print('\033[92m'+'\nEND OF PROGRAM'+'\033[0m')
