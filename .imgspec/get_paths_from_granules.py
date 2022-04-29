#! /usr/bin/python

import argparse
import glob
import os
import sys
import tarfile


def parse_args():
    parser = argparse.ArgumentParser()
    products = ["rfl"]
    formats = ["envi"]
    parser.add_argument("-p", "--product",
                        help=("Choose one of the following product types: " + ", ".join(products)))
    parser.add_argument("-f", "--format",
                        help=("Choose one of the following formats: " + ", ".join(formats)))
    args = parser.parse_args()

    if args.product:
        if args.product not in products:
            print("ERROR: Product \"%s\" is not a valid product choice." % args.product)
            sys.exit(1)
    if args.format:
        if args.format not in formats:
            print("ERROR: Format \"%s\" is not a valid format choice." % args.format)
            sys.exit(1)
    return args


def main():
    args = parse_args()

    # Unzip and untar granules
    input_dir = "input"
    granule_paths = glob.glob(os.path.join(input_dir, "*.tar.gz"))
    for g in granule_paths:
        tar_file = tarfile.open(g)
        tar_file.extractall(input_dir)
        tar_file.close()
        os.remove(g)

    # Get paths based on product type file matching
    paths = []
    if args.product == "rfl":
        # AVIRIS SDS uses *rfl*img and *corr*img for distributed files
        paths = glob.glob(os.path.join(input_dir, "*", "*rfl*img"))
        paths += glob.glob(os.path.join(input_dir, "*", "*corr*img"))
        # ISOFIT uses *rfl for processed reflectance files
        paths += glob.glob(os.path.join(input_dir, "*", "*rfl"))

    print(",".join(paths))


if __name__ == "__main__":
    main()
