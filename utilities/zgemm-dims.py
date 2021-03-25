#!/usr/bin/env python3

import os # to check file existence
import sys # for argc, argv
import re # regex
import lzma as xz
import numpy as np
import pandas as pd

# Define help text function
argv = sys.argv;
def help():
    print();
    print( argv[0], '- script to process Exciting-Plus ZGEMM/3M dimensions'   );
    print(  'Usage:', argv[0], '[input.xz]'                                   );
    print(  'where input.xz is an LZMA-compressed text file containing the'   );
    print(  '  output of the command'                                         );
    print( r'    grep -n " zgemm:" [LSF logfile]'                             );
    print(  '  with the LSF job log file from a response or c-RPA job run'    );
    print( r'  on Summit with a debug version (DEBUG=1) of Exciting-Plus.'    );
    print();

# Quick exit if no input file given
if len(argv) < 2:
    help();
    sys.exit( 'Error: No input file' );

# Quick exit if too many input files given
if len(argv) > 2:
    help();
    sys.exit( 'Error: Too many input files' );

# Quick exit if file does not exist
if not os.path.exists(argv[1]):
    help();
    sys.exit( 'Error: cannot read file', argv[1] );
else:
    
    # Define column names and columns to take
    cols     = [ 4,   7,   10  ];
    colnames = [ 'm', 'n', 'k' ];

    # Read file through LZMA decompressor
    ifname = argv[1];
    ifile = xz.open( ifname, 'rt' );

    # Read data line by line, filtering through regex
    # Adapted from https://stackoverflow.com/questions/66767377/filter-out-non-numeric-lines-from-multiple-columns-in-pandas/66768776
    pattern = r'\d+:\szgemm:\s+m\s=\s+(\d+)\s+n\s=\s+(\d+)\s+k\s=\s+(\d+)\s+$';
    rx = re.compile(pattern);
    data = pd.DataFrame(( m.groups() for line in ifile
                          for m in ( rx.match(line), ) if m ),
                          columns=colnames ).astype('int64');
    ifile.close();

    # Filter outliers
    # Adapted from https://nextjournal.com/schmudde/how-to-remove-outliers-in-data
    x = 0.00001 # percentile to remove
    datam = data['m'];
    datan = data['n'];
    datak = data['k'];
    m_in = datam.between(datam.quantile(x), datam.quantile(1-x));
    n_in = datan.between(datan.quantile(x), datan.quantile(1-x));
    k_in = datak.between(datak.quantile(x), datak.quantile(1-x));
    datam_clean = datam[m_in]
    datan_clean = datan[n_in]
    datak_clean = datak[k_in]
    outliers = np.count_nonzero(m_in == False) \
               + np.count_nonzero(n_in == False) \
               + np.count_nonzero(k_in == False)
    print( 'Removed', str(outliers), 'data points' );
    
    # Find minimum and maximum points
    mdims = datam_clean.to_numpy();
    mmin = np.amin(mdims);
    mmax = np.amax(mdims);
    ndims = datan_clean.to_numpy();
    nmin = np.amin(ndims);
    nmax = np.amax(ndims);
    kdims = datak_clean.to_numpy();
    kmin = np.amin(kdims);
    kmax = np.amax(kdims);

    # Display output
    print( re.sub( ifname, '.xz', '' ), ':' );
    print( 'M =', mmin, '-', mmax );
    print( 'N =', nmin, '-', nmax );
    print( 'K =', kmin, '-', kmax );

    sys.exit(0);
