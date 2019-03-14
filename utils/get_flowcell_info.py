#!/usr/bin/env python

import h5py, argparse

# Gather our code in a main() function
def main(args):
    # Read fast5 file and print flowcell, kit, and platform information
    print(f"reading fast5 file: {args.fast5}")

    f = h5py.File(args.fast5, 'r')
    ti_attrs = f["UniqueGlobalKey"]["tracking_id"].attrs
    print(f"Platform  {ti_attrs['device_type']}")
    print(f"Flowcell ID  {ti_attrs['flow_cell_id']}")

    ct_attrs = f["UniqueGlobalKey"]["context_tags"].attrs
    print(f"Sequencing kit  {ct_attrs['sequencing_kit']}")
    print(f"Flowcell Type  {ct_attrs['flowcell_type']}")


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Read fast5 file and print flowcell, kit, and platform information.")

    # Script parameters
    parser.add_argument(
                    "fast5",
                    help = "pass fast5 file to the program",
                    metavar = "fast5")
    args = parser.parse_args()
  
    main(args)
