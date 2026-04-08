#!/usr/bin/env python3
"""
mask_cox1_artifact.py
=====================
Replace a known homopolymer artifact in an assembled COX1 minichromosome
with Ns, producing a GenBank-submittable sequence.

Background
----------
The assembled contig contains a run of Cs (CCCCCTCACCCCCTCCCCC) between two
well-supported flanking sequences. All 20 Eustrongylides Sanger reference
sequences either end at or show conflicting indels through this region,
consistent with an ONT homopolymer overcalling artifact that neither
fast- nor SUP-basecalling resolves. Replacing the uncertain run with Ns is
the scientifically conservative option: it reports that sequence is present
but its exact composition is unknown.

The flanking anchors are searched on the gene's strand (minus strand for
this assembly). Everything between the left anchor and right anchor is
replaced with N_COUNT Ns. The anchors themselves are preserved.

Usage
-----
    python3 mask_cox1_artifact.py \\
        -i assembly.fasta \\
        -o assembly_masked.fasta \\
        [-n 10]                       # number of Ns (default: 10)
        [-l TATGGAAA]                 # left anchor (default shown)
        [-r TTAACATGTTGAACT]          # right anchor (default shown)
        [--contig contig_1]           # contig to edit (default: first)
        [--strand minus]              # strand the gene is on (default: minus)

Output
------
- assembly_masked.fasta : original assembly with the artifact replaced by Ns
- Prints a summary showing the replaced region and coordinates
"""

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def find_on_strand(seq, anchor, strand):
    """Return 0-based position of anchor on the given strand."""
    if strand == "plus":
        pos = seq.find(anchor)
    else:
        rc_anchor = str(Seq(anchor).reverse_complement())
        pos_rc = seq.find(rc_anchor)
        if pos_rc < 0:
            return -1
        # Convert RC position back to a plus-strand position range
        # (pos_rc is where the RC of the anchor sits on the plus strand)
        # We return the plus-strand end of the anchor (right side), which
        # corresponds to the left anchor's position on the minus strand read
        # left-to-right.  Caller uses this to define the replacement window.
        return pos_rc
    return pos


def mask_assembly(in_fasta, out_fasta, left_anchor, right_anchor,
                  n_count, target_contig, strand):
    records = list(SeqIO.parse(in_fasta, "fasta"))
    if not records:
        sys.exit("ERROR: No sequences in input FASTA.")

    # Pick target contig
    if target_contig:
        matches = [r for r in records if r.id == target_contig]
        if not matches:
            ids = ", ".join(r.id for r in records)
            sys.exit(f"ERROR: Contig '{target_contig}' not found. Available: {ids}")
        rec = matches[0]
        others = [r for r in records if r.id != target_contig]
    else:
        rec = records[0]
        others = records[1:]

    seq = str(rec.seq)
    n = len(seq)

    # -------------------------------------------------------------------------
    # Find anchors.  For minus-strand gene: the gene is read right-to-left on
    # the plus strand.  We search for each anchor's reverse complement on the
    # plus strand, then define the replacement window in plus-strand space.
    # -------------------------------------------------------------------------
    if strand == "minus":
        la_rc = str(Seq(left_anchor).reverse_complement())   # left anchor RC
        ra_rc = str(Seq(right_anchor).reverse_complement())  # right anchor RC

        # On plus strand the order is: ...ra_rc... <gap> ...la_rc...
        # because on minus strand reading direction is right→left on plus strand
        la_pos = seq.find(la_rc)
        ra_pos = seq.find(ra_rc)

        if la_pos < 0:
            sys.exit(f"ERROR: Left anchor RC '{la_rc}' not found in {rec.id}.")
        if ra_pos < 0:
            sys.exit(f"ERROR: Right anchor RC '{ra_rc}' not found in {rec.id}.")

        # On minus strand: left_anchor appears AFTER right_anchor in plus coords
        # Window to replace = between end of ra_rc and start of la_rc
        replace_start = ra_pos + len(ra_rc)   # first base after right anchor RC
        replace_end   = la_pos                # last base before left anchor RC

        if replace_start >= replace_end:
            # Might be the other orientation; try swapping
            replace_start = la_pos + len(la_rc)
            replace_end   = ra_pos
            if replace_start >= replace_end:
                sys.exit(
                    f"ERROR: Anchors found but window is empty or inverted.\n"
                    f"  left_anchor RC '{la_rc}' at {la_pos}\n"
                    f"  right_anchor RC '{ra_rc}' at {ra_pos}\n"
                    f"  Computed window: {replace_start}–{replace_end}"
                )

    else:  # plus strand
        la_pos = seq.find(left_anchor)
        ra_pos = seq.find(right_anchor)
        if la_pos < 0:
            sys.exit(f"ERROR: Left anchor '{left_anchor}' not found in {rec.id}.")
        if ra_pos < 0:
            sys.exit(f"ERROR: Right anchor '{right_anchor}' not found in {rec.id}.")
        replace_start = la_pos + len(left_anchor)
        replace_end   = ra_pos

    if replace_start >= replace_end:
        sys.exit(
            f"ERROR: Replacement window is empty (start={replace_start}, end={replace_end}).\n"
            "Check that the anchors flank the artifact in the correct order."
        )

    original = seq[replace_start:replace_end]
    ns = "N" * n_count
    new_seq = seq[:replace_start] + ns + seq[replace_end:]

    # -------------------------------------------------------------------------
    # Report
    # -------------------------------------------------------------------------
    print(f"Contig : {rec.id} ({n} bp, gene on {strand} strand)")
    print(f"Left anchor  : {left_anchor!r}")
    print(f"Right anchor : {right_anchor!r}")
    print(f"Replaced (plus-strand coords, 1-based): {replace_start+1}–{replace_end}")
    print(f"  Original sequence ({len(original)} bp): {original!r}")
    print(f"  Replaced with    ({n_count} bp): {'N'*n_count!r}")
    print(f"New contig length: {len(new_seq)} bp")

    # -------------------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------------------
    new_rec = SeqRecord(
        Seq(new_seq),
        id=rec.id,
        name=rec.name,
        description=rec.description + f" [artifact at {replace_start+1}-{replace_end} replaced with {n_count}xN]",
    )
    out_records = [new_rec] + others
    with open(out_fasta, "w") as fh:
        SeqIO.write(out_records, fh, "fasta")
    print(f"\nOutput written to: {out_fasta}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input",   required=True,  help="Input assembly FASTA")
    p.add_argument("-o", "--output",  required=True,  help="Output masked FASTA")
    p.add_argument("-n", "--n-count", type=int, default=10,
                   help="Number of Ns to insert (default: 10)")
    p.add_argument("-l", "--left-anchor",  default="TATGGAAA",
                   help="Sequence immediately LEFT of artifact on gene strand "
                        "(default: TATGGAAA)")
    p.add_argument("-r", "--right-anchor", default="TTAACATGTTGAACT",
                   help="Sequence immediately RIGHT of artifact on gene strand "
                        "(default: TTAACATGTTGAACT)")
    p.add_argument("--contig", default=None,
                   help="Contig ID to edit (default: first contig)")
    p.add_argument("--strand", choices=["plus", "minus"], default="minus",
                   help="Strand the gene is on (default: minus)")
    args = p.parse_args()

    mask_assembly(
        in_fasta=args.input,
        out_fasta=args.output,
        left_anchor=args.left_anchor,
        right_anchor=args.right_anchor,
        n_count=args.n_count,
        target_contig=args.contig,
        strand=args.strand,
    )


if __name__ == "__main__":
    main()
