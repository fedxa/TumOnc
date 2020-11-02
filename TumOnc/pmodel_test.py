import TumOnc as to
import TumOnc.pmodel as pm
import argparse
import scipy.stats as st


def main(fn, minnmut=100):
    bn = to.binned_ns_load(fn)

    pm.plot_pmodel_test(bn, minnmut)


if __name__ == "__main__":
    # execute only if run as a script
    parser = argparse.ArgumentParser(description='Run the analysis')
    parser.add_argument("--pmodel",
                        help="Uses previously saved pmodel instead of calculation",
                        default=False)
    parser.add_argument("--minnmut",
                        help="Minimal muation count for averaging window",
                        default=100)

    args = parser.parse_args()

    if args.pmodel:
        print('Loading pmodel from %s instead of estimating' % args.pmodel)
        main(args.pmodel,int(args.minnmut))
    else:
        print("specify the pmodel!")
        exit()
