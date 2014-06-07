#! /usr/bin/env python

from __future__ import division, print_function
import sys
from numpy import array, append, linspace, log
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from collect import collect
from scipy.stats.kde import gaussian_kde

def main():
    """\
    Script for producing a plot of the a probability density function from 
    the CI coefficients of a RAPTOR output (evb.out).  This mimics the concept
    proposed in the MS-EVB3 paper. 

    CHANGELOG
    
    9-20-2013 DWS v1.0
    Initial build.  Currently the code outputs data to screen or plot it.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-b', '--bins', help='Number of bins for generating plot.',
                        default=200)
    parser.add_argument('-f', '--freeenergy', help='Use the probability density function '
                        'to generate a free energy profile from the largest MS-EVB amplitude.',
                        action='store_true', default=False)
    parser.add_argument('-fd', '--freeenergydiff', help='Use the probability density function '
                        'difference between the largest and second largest MS-EVB amplitude '
                        'to generate a free energy profile.',
                        action='store_true', default=False)
    parser.add_argument('files', help='The files that you want converted.'
                        , nargs='+')
    args = parser.parse_args()

    # Store command line arguments as convenient variables
    fh = args.files[0]
    nbins = args.bins
    genfe = args.freeenergy
    genfed = args.freeenergydiff

    # Collect EVB simulation data 
    evb_data = collect(fh)

    # Get the largest (first) and second largest (second) CI vector
    # from an MS-EVB simulation 
    tmp = [i[1] for i in evb_data['CI_VECTOR']]

    first = array([])
    second =array([])
    for i in tmp:
        l = max(i)
        l2 = second_largest(i)
        # Note that we store the squared CI coefficient
        first = append(first, l*l)
        second = append(second, l2*l2)

    # Generate probability densities.  Here we use the kernel density estimator
    # to give a smooth PDF.
    pdf1 = gaussian_kde(first)
    #bin_pdf1 = linspace( min(first), max(first), nbins)
    bin_pdf1 = linspace( 0, 1, nbins )
    pdf2 = gaussian_kde(second)
    #bin_pdf2 = linspace( min(second), max(second), nbins)
    bin_pdf2 = linspace( 0, 1, nbins )

    if genfed:
        pdf3 = gaussian_kde(first-second)
        pdf4 = gaussian_kde(second-first)
        bin_pdf3 = linspace( -1, 1, 2*nbins )

    # If requested, generate the free energy corresponding to the probability
    # density
    if genfe or genfed:
        free_energy = array([])
        kb = 1.3806488e-23 # J/K
        # 1 kcal/mol = 6.9477e-21 J
        kb_kcal = kb / 6.9477e-21
        if genfe:
            for i in bin_pdf1:
                val = -kb_kcal*300.0*log(pdf1(i))
                free_energy = append(free_energy, val)
        else:
            for i in bin_pdf3:
                if i < -0.04:
                    val = -kb_kcal*300.0*log(pdf4(i))
                elif i > 0.04:
                    val = -kb_kcal*300.0*log(pdf3(i))
                else:
                    val = 0.0
                free_energy = append(free_energy, val)
        # Shift the zero point on the graph.
        minfe = min(free_energy)
        free_energy += -minfe

    # Plot the data
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif', 'serif': ['Times'], 'size': 30})
    fig = plt.figure()
    sub = fig.add_subplot(111)
    if genfe:
        sub.plot(bin_pdf1, free_energy, 'b', linewidth=3)
    elif genfed:
        sub.plot(bin_pdf3, free_energy, 'b', linewidth=3)
    else:
        sub.plot(bin_pdf1, pdf1(bin_pdf1), 'b', linewidth=3)
        sub.plot(bin_pdf2, pdf2(bin_pdf2), 'r', linewidth=3)

    # Adjust graph scale
    fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

    # Title, labels
    if genfe:
        sub.set_xlabel(r'c$_1^2$')
        sub.set_ylabel(r'Free Energy (kcal/mol)')
    elif genfed:
        sub.set_xlabel(r'c$_1^2$-c$_2^2$')
        sub.set_ylabel(r'Free Energy (kcal/mol)')
    else:
        sub.set_xlabel(r'c$_i^2$')
        sub.set_ylabel(r'Probability Density')

    # Axis limits
    if genfe:
        sub.set_xlim([0.35,0.85])
        sub.set_ylim([min(free_energy),3])
    elif genfed:
        sub.set_xlim([-0.80,0.80])
        sub.set_ylim([min(free_energy),3])
    else:
        sub.set_xlim([0,0.85])
        sub.set_ylim([0,8])

    # Minor tick marks
    if genfed:
        xminor = MultipleLocator(0.05)
    else:
        xminor = MultipleLocator(0.01)
    yminor = MultipleLocator(0.1)
    sub.xaxis.set_minor_locator(xminor)
    sub.yaxis.set_minor_locator(yminor)

    plt.show()

def second_largest(numbers):
    '''Function for determining the second largest element in an array.'''
    m1, m2 = None, None
    for x in numbers:
        if x >= m1:
            m1, m2 = x, m1
        elif x > m2:
            m2 = x
    return m2
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
