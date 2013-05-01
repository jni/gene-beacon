import argparse
import numpy as np
import progressbar
import Bio
import Bio.SeqIO
from Bio.pairwise2 import align

progress_widgets = [progressbar.Percentage(), ' ', progressbar.Counter(), ' ',
                    progressbar.ETA(), ' ', progressbar.Bar(marker='=')]


def get_sequence_string(seq):
    """Extract just the sequence contained in a SeqIO record as a string.

    Parameters
    ----------
    seq : Bio.SeqRecord or string
        The sequence to be converted.

    Returns
    -------
    seqstr : string
        A string containing the sequence.
    """
    if type(seq) == Bio.SeqRecord:
        seqstr = seq.seq.tostring()
    elif type(seq) == Bio.Seq.Seq:
        seqstr = seq.tostring()
    else:
        seqstr = seq
    return seqstr


def match_score(seq1, seq2):
    """Match the two input sequences using local alignment.

    Parameters
    ----------
    seq1 : string or Bio.SeqRecord
    seq2 : string or Bio.SeqRecord
        The two sequences to be matched.

    Returns
    -------
    score : float
        The score of the alignment.
    """

    seq1 = get_sequence_string(seq1)
    seq2 = get_sequence_string(seq2)
    score = align.localxx(seq1, seq2)[0][2]
    return score


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Classify sequences as matching one or'
        ' another probe.')
    parser.add_argument('-p', '--probe-sequences', metavar='FILENAME',
        help='The file containing the probe sequences.',
        required=True)
    parser.add_argument('-s', '--test-sequences', metavar='FILENAME',
        help='The file containing the test sequences.',
        required=True)
    parser.add_argument('-f', '--probe-file-format', default='fasta',
        help='The format of the sequences in the probe file. (default: fasta)')
    parser.add_argument('-F', '--test-file-format', default='fastq',
        help='The format of the sequences in the test file. (default: fastq)')
    parser.add_argument('-P', '--progress', default=False, action='store_true',
        help='Report progress with an ASCII progress bar.')

    args = parser.parse_args()
    with open(args.probe_sequences) as probe_file:
        probes = list(Bio.SeqIO.parse(probe_file, args.probe_file_format))
    with open(args.test_sequences) as test_file:
        test = list(Bio.SeqIO.parse(test_file, args.test_file_format))

    out_lists = [[] for i in range(len(probes))]
    if args.progress:
        pbar = progressbar.ProgressBar(widgets=progress_widgets)
        test = pbar(test)
    for record in test:
        scores = [match_score(record.seq, pr) for pr in probes]
        out_lists[np.argmax(scores)].append(record)

    for i, group in enumerate(out_lists):
        fn = args.test_sequences + '.%d'%i
        with open(fn, 'w') as fout:
            Bio.SeqIO.write(group, fout, args.test_file_format)
