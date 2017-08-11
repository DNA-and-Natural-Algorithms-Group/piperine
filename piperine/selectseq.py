
# E.g. run like this:
# python selectseq.py oscillator_scores_bmax10.csv
# python selectseq.py oscillator_scores_bmax2000.csv
# python selectseq.py oscillator_scores_bmax10.csv oscillator_scores_bmax2000.csv

import sys, csv
import numpy

if __name__ == "__main__":

    csvnames = sys.argv[1:]

    columns_named=False
    scores = []
    for csvname in sys.argv[1:] :
        print("Reading scores from "+csvname+" .")

        with open(csvname, 'rb') as csvfile:
            score_reader = csv.reader(csvfile)
            for row in score_reader:
                print(', '.join(row))
                if 'Index' in row[0] and columns_named:
                    continue
                else:
                    scores.append(row)
                    columns_named=True
            csvfile.close()

    numseqs = len(scores)-1
    numscores = len(scores[1])-1
    print("Found information for " + str(numseqs) + " sequences with " + str(numscores) + " score types each.")

    columns = list(zip(*scores))
    ranks = []
    fractions = []
    percents = []
    for col in columns:
        # print col
        if 'Index' in col[0] or 'Defect' in col[0] or 'Toehold Avg' in col[0] or 'Range of toehold' in col[0]:
            continue
        if 'SSU' in col[0] or 'SSTU' in col[0]:    # for these scores, higher is better
            col = [-float(x) for x in col[1:]]
        else:
            col = [float(x) for x in col[1:]]
        array = numpy.array(col)
        temp = array.argsort()
        colranks = numpy.empty(len(array), int)
        colranks[temp] = numpy.arange(len(array))
        ranks.append(colranks)                     # low rank is better
        fractions.append((array - array.min())/abs(array.min() + (array.min()==0) ))
        percents.append((array - array.min())/(array.max() - array.min()))
    temp=ranks
    ranks=list(zip(*temp))
    temp=fractions
    fractions=list(zip(*temp))
    temp=percents
    percents=list(zip(*temp))

    print("\nRank array:")
    print("                         ")
    for title in scores[0]:
        if 'Index' in title or 'Defect' in title or 'Toehold Avg' in title or 'Range of toehold' in title:
            continue
        print("{:>6s}".format(title[0:6]))
    print()
    i=0
    for r in ranks:
        print("design {:2d}: {:6d} = sum [".format(i,sum(r)))
        for v in r:
            print("{:6d}".format(v))
        print("]")
        i=i+1

    print("\nFractional excess array:")
    print("                         ")
    for title in scores[0]:
        if 'Index' in title or 'Defect' in title or 'Toehold Avg' in title or 'Range of toehold' in title:
            continue
        print("{:>6s}".format(title[0:6]))
    print()
    i=0
    for f in fractions:
        print("design {:2d}: {:6.2f} = sum [".format(i,sum(f)))
        for v in f:
            print("{:6.2f}".format(v))
        print("]")
        i=i+1

    print("\nPercent badness (best to worst) array:")
    print("                         ")
    for title in scores[0]:
        if 'Index' in title or 'Defect' in title or 'Toehold Avg' in title or 'Range of toehold' in title:
            continue
        print("{:>6s}".format(title[0:6]))
    print()
    i=0
    for p in percents:
        print("design {:2d}: {:6.2f} = sum [".format(i,100*sum(p)))
        for v in p:
            print("{:6.2f}".format(100*v))
        print("]")
        i=i+1

    print(" ")
    worst_rank = 0
    while 1:
        ok_seqs = [i for i in range(len(ranks)) if max(ranks[i])<=worst_rank]
        if len(ok_seqs)==0:
            worst_rank=worst_rank+1
            continue
        else:
            break

    # scores used:
    # TSI avg, TSI max, TO avg, TO max, BM, Largest Match, SSU Min, SSU Avg, SSTU Min, SSTU Avg, Max Bad Nt %,  Mean Bad Nt %, WSI-Intra, WSI-Inter, WSI-Intra-1, WSI-Inter-1, Verboten, WSI
    weights = [5,   20,     10,     30,  2,             3,      30,      10,       50,       20,           10,              5,         6,         4,           5,           3,        2,  8]

    print("Indices of sequences with best worst rank of " + str(worst_rank) + ": " + str(ok_seqs))
    print("  Sum of all ranks, for these sequences:      " + str([sum(ranks[i]) for i in ok_seqs]))
    print("  Sum of weighted ranks, for these sequences: " + str([sum(numpy.array(ranks[i])*weights/100.0) for i in ok_seqs]))
    print("  Sum of fractional excess over best score:   " + str([sum(fractions[i]) for i in ok_seqs]))
    print("  Sum of weighted fractional excess:          " + str([sum(numpy.array(fractions[i])*weights/100.0) for i in ok_seqs]))
    print("  Sum of percent badness scores:              " + str([100*sum(percents[i]) for i in ok_seqs]))
    print("  Sum of weighted percent badness scores:     " + str([sum(numpy.array(percents[i])*weights) for i in ok_seqs]))
    temp = [sum(r) for r in ranks]
    print("Best sum-of-ranks:                   {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), numpy.argmin(temp), max(temp), numpy.argmax(temp) ))
    temp = [sum(numpy.array(r)*weights/100.0) for r in ranks]
    print("Best sum-of-weighted-ranks:          {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), numpy.argmin(temp), max(temp), numpy.argmax(temp) ))
    temp = [sum(f) for f in fractions]
    print("Best fractional excess sum:          {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), numpy.argmin(temp), max(temp), numpy.argmax(temp) ))
    temp = [sum(numpy.array(f)*weights/100.0) for f in fractions]
    print("Best weighted fractional excess sum: {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), numpy.argmin(temp), max(temp), numpy.argmax(temp) ))
    temp = [100*sum(p) for p in percents]
    print("Best percent badness sum:            {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), numpy.argmin(temp), max(temp), numpy.argmax(temp) ))
    temp = [sum(numpy.array(p)*weights) for p in percents]
    print("Best weighted percent badness sum:   {:6.2f} by [{:d}]      and the worst: {:6.2f} by [{:d}]".format( min(temp), numpy.argmin(temp), max(temp), numpy.argmax(temp) ))
    print()
