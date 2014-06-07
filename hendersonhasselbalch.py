#! /usr/bin/env python

from __future__ import division, print_function
import sys
from numpy import arange, append, array, sort
import matplotlib.pyplot as plt

def main():
    """\
    Script for evaluating the Henderson-Hasselbalch equation.
    """

    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    # Carbonic acid pKa = 3.45 (Ka = 2.5 x 10^{-4})
    # Bicarbonate pKa = 10.329 (Ka = 4.69 x 10^{-11})
    pka1 = 3.45
    pka2 = 10.329

    ka1 = 3.54e-4
    ka2 = 4.69e-11

    # Mass balance and acid equilibria expressions allow us to write:
    # H2CO3 fraction = [H3O+]^2 / ( [H3O+]^2 + [H3O+]*Ka1 + Ka1*Ka2 )
    # HCO3 fraction = [H3O+]*Ka1 / ( [H3O+]^2 + [H3O+]*Ka1 + Ka1*Ka2 )
    # CO3 fraction = Ka1*Ka2 / ( [H3O+]^2 + [H3O+]*Ka1 + Ka1*Ka2 )
 
    # pH values
    phvals = arange(0.0,14.1,0.1)
    phvals = append(phvals, pka1)
    phvals = append(phvals, pka2)
    phvals = sort(phvals)

    # Fractions of species
    h2co3 = array([])
    hco3 = array([])
    co3 = array([])

    # pH values
    pp = ( 1.0, 2.0, 3.0, pka1, 4.0, 5.0,
           6.0, 7.0, 8.0, 9.0, 10.0, pka2,
           11.0, 12.0, 13.0, 14.0 
         )

    fmt = '{0:6.3f} {1:8.6f} {2:8.6f} {3:8.6f}'

    print()
    print('   pH   H2CO3    HCO3^-  CO3^{2-}')
    for ph in phvals:
        h3o = 10**(-ph)
        tmp1 = h3o*h3o / ( h3o*h3o + h3o*ka1 + ka1*ka2 )
        tmp2 = h3o*ka1 / ( h3o*h3o + h3o*ka1 + ka1*ka2 )
        tmp3 = ka1*ka2 / ( h3o*h3o + h3o*ka1 + ka1*ka2 )
        h2co3 = append(h2co3, tmp1)
        hco3 = append(hco3, tmp2)
        co3 = append(co3, tmp3)
        if ph in pp:
            print(fmt.format(ph, tmp1, tmp2, tmp3))
    print()

    # Make a plot
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif', 'serif': ['Times'], 'size': 30})

    fig = plt.figure()
    sub = fig.add_subplot(111)
    sub.plot(phvals, h2co3, 'b', linewidth=3)
    sub.plot(phvals, hco3, 'g', linewidth=3)
    sub.plot(phvals, co3, 'r', linewidth=3)

    # Draw vertical lines
    sub.plot([pka1, pka1], [0, 1], 'k--', lw=2)
    sub.plot([pka2, pka2], [0, 1], 'k--', lw=2)

    fig.subplots_adjust(left=0.12, right=0.90, bottom=0.1, top=0.9)

    # Title, labels
    sub.set_xlabel(r'pH')
    sub.set_ylabel(r'Fraction of species')

    # Limits
    sub.set_ylim([0,1.1])

    # Tick mark labels
    sub.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14])
    sub.set_yticks([0,0.2,0.4,0.6,0.8,1.0])

    # Text
    sub.text(0.5, 1.02, r'H$_2$CO$_3$')
    sub.text(3.2, 1.02, r'pK$_{a1}$')
    sub.text(6.5, 1.02, r'HCO$_3^-$')
    sub.text(10.05, 1.02, r'pK$_{a2}$')
    sub.text(12.5, 1.02, r'CO$_3^{2-}$')

    plt.show()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
