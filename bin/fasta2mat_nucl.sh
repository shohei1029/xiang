#!/bin/env bash
set -euo pipefail

# nucleotide only

# 2017.9.5 by Shohei Nagata
# Multi FASTAを読み込ませると，sim, matを出力する.

# 9.12
# 既に作成されていたblastp_makesim_args_v3.shを組み込む形で。。。。あったんかいェ....。

# 10.19
# xiangライブラリ内に移行。合わせてblast2sim.pyコマンドがpathに入ってるかの判定でエラー表示を出すように。

EVALUE=1e-5
#SCRIPT_BLAST_TO_SIM=blast2sim.py

FASTA_FILE=$1

if [[ ! -f ${FASTA_FILE} ]]; then
    echo "can't find ${FASTA_FILE}, exitting.."
    exit 1
fi

if [ ! `which blast2sim.py` ];then
    echo "can't find blast2sim.py in your path, exitting.."
    exit 1
fi

NAME_CORE=`basename ${FASTA_FILE} .fasta`
OUTDIR=out_${NAME_CORE}
mkdir -p ${OUTDIR}


#########
# BLAST #
#########
echo "making BLAST Database.."
makeblastdb -in ${FASTA_FILE} -out ${OUTDIR}/${NAME_CORE} -dbtype nucl -hash_index -parse_seqids -max_file_sz 2GB

BLAST_OUT_FILE_NAME=blastn_bestHSP_1e-5_${NAME_CORE}.txt

echo "BLAST.."
blastn -query ${FASTA_FILE} -db ${OUTDIR}/${NAME_CORE} -max_hsps 1 -num_threads 16 -outfmt 6 -evalue ${EVALUE} -max_target_seqs 99999999 -out ${OUTDIR}/${BLAST_OUT_FILE_NAME}

###################
# create sim file #
###################

#-max_hsps 1 オプションと，-outfmt 6により不要に
#grep -v "#" ${OUTDIR}/${BLAST_OUT_FILE_NAME} | LC_ALL=C sort -k 1,2 -u > ${OUTDIR}/mbs_${BLAST_OUT_FILE_NAME}

echo "creating sim files.."
blast2sim.py -i ${OUTDIR}/${BLAST_OUT_FILE_NAME} > ${OUTDIR}/sim_${BLAST_OUT_FILE_NAME}

###
jobline "$0 ${NAME_CORE} done"
