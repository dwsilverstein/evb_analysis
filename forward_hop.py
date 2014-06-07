#! /usr/bin/env python

from __future__ import division, print_function
import sys
from numpy import array
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from collect import collect

def main():
    """\
    Script for producing a plot of the forward hopping of protons from 
    a RAPTOR output (evb.out).  This is only tested for the case of 
    a dissolved proton (e.g. H3O in a box of nH2O). 

    CHANGELOG
    
    9-20-2013 DWS v1.0
    Initial build.  Currently the code outputs data to screen or plot it.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('files', help='The files that you want converted.'
                        , nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files[0]

    # Collect EVB simulation data 
    evb_data = collect(fh)

    # Get the timestep and rxncenters from an MS-EVB simulation 
    ts = evb_data['TIMESTEP']
    rxncenter = [i[1] for i in evb_data['RXNCENTER']]

    hopfxn = eval_hop_function(ts, rxncenter)

    # Plot the data
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif', 'serif': ['Times'], 'size': 30})
    fig = plt.figure()
    x = array(ts)/1000.0
    y = hopfxn
    sub = fig.add_subplot(111)
    sub.plot(x, y, 'b', linewidth=3)

    # Adjust graph scale
    fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

    # Title, labels
    sub.set_xlabel(r'time (ps)')
    sub.set_ylabel(r'Forward Hop')

    # Axis limits
    sub.set_xlim([0,x[-1]])
    #sub.set_xlim([1300,1600])
    sub.set_ylim([0,y[-1]])
    #sub.set_ylim([480,610]) 

    plt.show()

def eval_hop_function(ts, rxncenter):
    '''Evaluate the hop function from the MS-EVB3 paper.  This is defined:
       
       h(ts) = h(ts-1) + dh(ts)
       h(0) = 0
    
                (  0 , if no proton hop
       dh(ts) = {  1 , if proton hops to new acceptor
                ( -1 , if proton hops to the last donor

       Here, ts is the timestep, and ts-1 is the previous timestep.'''

    # Initialize the hop function
    hts = []
    # Set h(0) = 0
    hts.append(0)

    # Number of timesteps
    numts = len(ts)

    # Initialize the identity of the previous donor
    pdonor = [0]
    for i in range(1,numts):
        # Previous center
        pcenter = rxncenter[i-1]
        # Current center
        ccenter = rxncenter[i]
        # Evaluate the function
        if pcenter == ccenter:
            # No hop occurred
            val = 0
        else:
            # Hop occurred
            if ccenter == pdonor[-1]:
                # Here we checked that the current center is the previous donor.
                # If it is, we hopped backwards.
                val = -1
            else:
                # If the current center changed (ccenter != pcenter) and the
                # current center is not the previous donor (ccenter != pdonor[-1]),
                # we must have hopped forward.
                val = 1
                # Only change the last element of pdonor if a previous forward hop occurred.
                pdonor.append(pcenter) 
        hts.append(hts[i-1] + val)
    
    return hts
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
