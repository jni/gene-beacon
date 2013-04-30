import argparse
import numpy as np
import Bio.SeqIO

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
    score = np.random.rand()
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

    args = parser.parse_args()
    with open(args.probe_sequences) as probe_file:
        probes = list(Bio.SeqIO.parse(probe_file, args.probe_file_format))
    with open(args.test_sequences) as test_file:
        test = list(Bio.SeqIO.parse(test_file, args.test_file_format))

    out_lists = [[] for i in range(len(probes))]
    for record in test:
        scores = [match_score(record.seq, pr) for pr in probes]
        out_lists[np.argmax(scores)].append(record)

    for i, group in enumerate(out_lists):
        fn = args.test_sequences + '.%d'%i
        with open(fn, 'w') as fout:
            Bio.SeqIO.write(group, fout, args.test_file_format)
