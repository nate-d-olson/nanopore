import numpy
import pysam
import sys
import tqdm

def n50(values, fraction=0.5):
    if len(values) < 5:
        return numpy.nan
    values = values.copy()
    values.sort()
    
    cumsum = numpy.cumsum(values)
    sum_ = cumsum[-1]
    
    i = numpy.where(cumsum>=fraction*sum_)[0][0]
    
    return values[i]
    

def coverages(lengths, cutoffs):
    covs = []
    for cutoff in cutoffs:
        cur_lengths = lengths[lengths >= cutoff]
        covs.append(cur_lengths.sum() / 3.1e9)
    return covs
    

def plot_coverages(path, lengths):
    try:
        from biorpy import r
    except ImportError:
        print("biorpy isn't installed; figures won't be generated (install from https://github.com/nspies/biorpy)")
        return

    r.pdf(path)
    x = numpy.arange(0,0.5e6+1,1e4)
    y = coverages(lengths, x)

    r.plot(x/1e3, y, type="l", xlim=[0,575], lwd=2, cex=1.5,
           xlab="read length (kb)", ylab="coverage by reads > length", main=f"coverage by mapped reads")

    highlight_x = [0, 50e3, 100e3, 250e3, 500e3]
    highlight_y = coverages(lengths, highlight_x)
    r.points(numpy.array(highlight_x)/1e3, highlight_y, pch=20)
    r.text(numpy.array(highlight_x)/1e3, highlight_y, 
           ["{:.2f} ({:.1%})".format(i, i/highlight_y[0]) for i in highlight_y], pos=4)
    r.mtext(f"{lengths.sum()/1e6:,.1f}mb total")

    r.devoff()


def do_qc(path, pdf_path=None):
    inf = pysam.AlignmentFile(path)
    
    read_lengths = []
    unmapped_count = 0
    unmapped_bases = 0

    for count, read in tqdm.tqdm(enumerate(inf)):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_unmapped:
            unmapped_bases += read.query_length
            unmapped_count += 1
            continue

        read_lengths.append(read.query_alignment_length)
        # if count > 2000:
        #     print("ENDING EARLY:", str(read)[:1000])
        #     break
        
    read_lengths = numpy.array(read_lengths)

    print(f"Number of mapped reads: {len(read_lengths):,} (excludes supplementary and secondary alignments)")
    print(f"Number of unmapped reads: {unmapped_count:,}")
    print(f"Number of unmapped bases: {unmapped_bases:,}")

    print(f"N50: {n50(read_lengths):,}")
    cutoffs = [0, 10e3, 25e3, 50e3, 100e3, 250e3, 500e3]
    covs = coverages(read_lengths, cutoffs)
    for cutoff, cov in zip(cutoffs, covs):
        print(f" coverage by reads >= {int(cutoff):>10,}: {cov:.3f}x ({cov/covs[0]:6.1%})")
    
    print("Top read lengths:")
    read_lengths.sort()
    
    for length in read_lengths[:-11:-1]:
        print(f" {length:,}")
    
    if pdf_path:
        plot_coverages(pdf_path, read_lengths)


if __name__ == "__main__":
    do_qc(*sys.argv[1:])
