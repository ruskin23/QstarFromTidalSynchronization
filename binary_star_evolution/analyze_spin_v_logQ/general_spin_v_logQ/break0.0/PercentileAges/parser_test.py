import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m',
                    action='store',
                    dest='sampling_method',
                    help='select a sampling method for mcmc: uncorrelated or adaptive'
                    )

args=parser.parse_args()

if args.sampling_method:print(args.sampling_method)
