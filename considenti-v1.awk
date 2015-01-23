# Consensus Identity Index v1.0, 27 Dec 2014
# Calculating Consensus identity index to set of sequences
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
#
# v1.0, 27 Dec 2014 - initial release
#
# Required params: -v cons="consensus-align-seq.fas" - FASTA file with consensus sequense(s)
#                             subject-align-seq.fas  - FASTA file with subject sequense(s)
#
# Options: -v abovediag=1  - prints index values above diagonal in the top right corner
#          -v belowdiag=1  - prints index values below diagonal in the bottom left corner
#          -v similarity=1 - prints percent identity of the two compared sequences in opposed corner
#
# Usage: gawk -v cons="consensus-align-seq.fas" -f considenti-v1.awk subject-align-seq.fas > output.txt


BEGIN {
  # ignoring the letter case in sequences
  IGNORECASE = 1;
  # setting the default position to print index
  if ((!belowdiag)||(abovediag && belowdiag)) { abovediag = 1; belowdiag = 0 }

  # array preparation with query sequences
  while ((getline < cons) > 0) {
    if ($0 ~ /^>/) {
      seqHeader = $0;
      if (seq) {
        qSum++;
        qSeqA[qSum][prevSeqHeader][0] = "";
        gsub(/\n/, "", seq);
        qSeqLen = split(seq, qSeqA[qSum][prevSeqHeader], "");
      }
      prevSeqHeader = seqHeader;
      seq = "";
    }
    # appending multiline sequences
    if ($0 !~ /[^ACGTN-]/) { seq = seq $0 }
  }
  # for the last
  if (seq) {
    qSum++;
    qSeqA[qSum][prevSeqHeader][0] = "";
    gsub(/\n/, "", seq);
    qSeqLen = split(seq, qSeqA[qSum][prevSeqHeader], "");
    seq = "";
    prevSeqHeader = "";
  }
  close(cons);
}

# collecting monomers to array with subject sequences
/^>/ {
  seqHeader = $0;
  if (seq) {
    sSum++;
    sSeqA[sSum][prevSeqHeader][0] = "";
    gsub(/\n/, "", seq);
    sSeqLen = split(seq, sSeqA[sSum][prevSeqHeader], "");
  }
  prevSeqHeader = seqHeader;
  seq = "";
}
# appending multiline sequences
!/[^ACGTN-]/ { seq = seq $0 }

END {
  # for the last
  if (seq) {
    sSum++;
    sSeqA[sSum][prevSeqHeader][0] = "";
    gsub(/\n/, "", seq);
    sSeqLen = split(seq, sSeqA[sSum][prevSeqHeader], "");
    seq = "";
    prevSeqHeader = "";
  }
  # checking the length
  if (qSeqLen != sSeqLen) {
    print "Error: length of the query and the subject sequences not equal";
    exit;
  }

  print "Consensus Identity Index:\n";

  # Now we have two arrays with query - consensuses and subject sequenses
  # qSeqA[seq_number][">seq_header"][base_position]
  # sSeqA[seq_number][">seq_header"][base_position]

  # For each query we look over the subjects, and in the current position
  # glance down to base following sequences
  for (seqNumQuery=1; seqNumQuery<=qSum; seqNumQuery++) {
    for (seqHeaderQuery in qSeqA[seqNumQuery]) {
      for (seqNumSubject=1; seqNumSubject<=sSum; seqNumSubject++) {
        for (seqHeaderSubject in sSeqA[seqNumSubject]) {
          for (base=1; base<=qSeqLen; base++) {
            baseQuery = qSeqA[seqNumQuery][seqHeaderQuery][base];
            baseSubject = sSeqA[seqNumSubject][seqHeaderSubject][base];
            if (baseQuery != baseSubject) {
              # counting differences in sites with consensus and subject sequence
              diffQA[seqNumQuery][seqNumSubject]["baseDiff"]++;
            } else {
              # counting identical bases in consensus and subject sequence
              simQA[seqNumQuery][seqNumSubject]["baseSim"]++;
            }
            for (seqNumSubjectNext=seqNumSubject+1; seqNumSubjectNext<=sSum; seqNumSubjectNext++) {
              for (seqHeaderSubjectNext in sSeqA[seqNumSubjectNext]) {
                baseSubjectNext = sSeqA[seqNumSubjectNext][seqHeaderSubjectNext][base];
                if (baseQuery != baseSubjectNext) {
                  if (baseSubject == baseSubjectNext) {
                    # counting coinciding differences between subject and sequence below subject
                    coiDiffA[seqNumQuery][seqNumSubject][seqNumSubjectNext]["coiDiff"]++;
                  }
                }
                if (baseSubject == baseSubjectNext) {
                  # counting identical bases in subject and sequence below subject
                  simA[seqNumQuery][seqNumSubject][seqNumSubjectNext]["baseSim"]++;
                }
              }
            }
          }
        }
      }
    }
  }

  # printing separate tables for each query
  for (seqNumQuery=1; seqNumQuery<=qSum; seqNumQuery++) {
    for (seqHeaderQuery in qSeqA[seqNumQuery])
      printf("\n%s\n", seqHeaderQuery);
    for (seqNumSubject=1; seqNumSubject<=sSum; seqNumSubject++) {
      for (seqHeaderSubject in sSeqA[seqNumSubject])
        printf("\n%s", seqHeaderSubject);
      for (seqNumSubjectNext=1; seqNumSubjectNext<=sSum; seqNumSubjectNext++) {
        if (abovediag) {
          if (seqNumSubject < seqNumSubjectNext) {
            coiDiff = coiDiffA[seqNumQuery][seqNumSubject][seqNumSubjectNext]["coiDiff"];
            baseQDiff = diffQA[seqNumQuery][seqNumSubject]["baseDiff"];
            baseDiff = diffQA[seqNumQuery][seqNumSubjectNext]["baseDiff"];
            if (baseDiff) {
              printf("\t%.f", ((coiDiff*2)/(baseDiff+baseQDiff))*100);
            } else {
              printf("\t%d", 0);
            }
          }
          if (seqNumSubject == seqNumSubjectNext) printf("\t-");
          if (similarity) {
            if (seqNumSubject > seqNumSubjectNext) {
              baseSim = simA[seqNumQuery][seqNumSubjectNext][seqNumSubject]["baseSim"];
              printf("\t%.f", (baseSim/qSeqLen)*100);
            }
          } else if (seqNumSubject > seqNumSubjectNext) printf("\t%d", 0);
        }
        if (belowdiag) {
          if (seqNumSubject > seqNumSubjectNext) {
            coiDiff = coiDiffA[seqNumQuery][seqNumSubjectNext][seqNumSubject]["coiDiff"];
            baseQDiff = diffQA[seqNumQuery][seqNumSubjectNext]["baseDiff"];
            baseDiff = diffQA[seqNumQuery][seqNumSubject]["baseDiff"];
            if (baseDiff) {
              printf("\t%.f", ((coiDiff*2)/(baseDiff+baseQDiff))*100);
            } else {
              printf("\t%d", 0);
            }
          }
          if (seqNumSubject == seqNumSubjectNext) printf("\t-");
          if (similarity) {
            if (seqNumSubject < seqNumSubjectNext) {
              baseSim = simA[seqNumQuery][seqNumSubject][seqNumSubjectNext]["baseSim"];
              printf("\t%.f", (baseSim/qSeqLen)*100);
            }
          } else if (seqNumSubject < seqNumSubjectNext) printf("\t%d", 0);
        }
      }
    }
    printf("\n");
  }
}
