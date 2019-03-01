#!/ussr/bin/sh
# Generate tar balls with individual multi-sequence fast5 files consistent with copy_to_server.sh

# Run script from within the directory with fast5 files to tar
# Usage: ./tar_multiseq_fast5.sh RUN

## Pass Run from nanopore database as first argument
RUN=${1}

for fast5 in *fast5; do
	idx=${fast5##*_}
	idx=${idx%.*5}
	TARNAME=${RUN}_${idx}.tar
	tar -cvf ${TARNAME} ${fast5}
done
