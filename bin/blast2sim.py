#!/usr/bin/env python3

import re
import time
import sys
from collections import defaultdict
import argparse
import logging

import numpy as np

logger = logging.getLogger(__name__)
handler = logging.StreamHandler() #sys.stderr
handler.setLevel(logging.INFO)
logger.setLevel(logging.INFO)
logger.addHandler(handler)

# created by Shohei Nagata
# blastの結果ファイル (format: 6, tab-separated)を読み込ませる。予めmax bitscoreのみにしておく。
# sim (ID1 ID2 sim_score) を出力する.

# v3との変更点(change log)
# outpupt filename の変更
#  -> edge絞りの情報削除
# default で bothfmt
# --nomultihit も default で True (bit scoreの最大値のみ記録されているファイルを読み込ませる)

# * 類似度スコアでのcutoffを"未満"から"以下"へ変更 ***

# v4からの変更点(change log)
# scps用 sim 形式のみの出力へ (cytoscapeもその形式で読み込むし) <ID1 ID2 sim>
# Introduction of logging and progressbar.
# defaultで，sys.stdinから読み込んでsys.stdoutへ出力 (unix!)

# 16.3.1 BLAST idendityでのcutoffオプションの追加 #unix commandでやったほうが高速です

# v5naive
# progressbarがpipでなぜか入らなくなっているので，コメンアウト

# v6
# ついに大規模リニューアル！
# 1年生の時に書いたスクリプトの負の遺産を一掃。
# クラスとか使っちゃう♪

#16spp/scripts内にあるやつより新しい。

# 17.9.21
# 引数系修正
# まだ外部から呼び出すことにはうまく対応できてないかも e.g. logging系の呼び出し場所，引数系



###################
# Argument Parser #
###################
parser = argparse.ArgumentParser()
#parser.add_argument("-m","--nomultihit",type=bool,default=True,choices=[True,False],help="If there are NO multiple hits in Blast output file -> True .(分割ヒット対策をしなくていいんだったらTrue, 分割ヒットしててbit scoreの最大値を選んでとってくる必要があるんだったらFalse (Defauls)")
parser.add_argument("-i","--input_file",type=str,default=None,help="Input Blast output file (tabular).")
parser.add_argument("-o","--output_file",type=str,help="Output file")
#parser.add_argument("-k","--edge_num",type=int,default=0,help="Ideal number of edges. Conserve alledges -> 0（エッジ数絞らなくていいなら0を指定")
parser.add_argument("-w","--cutoff_edge_weight",type=float,default=0.0,help="Cutoff by edge weight (i.e. sim(x,y)).(<=),以下")
parser.add_argument("-s","--cutoff_identity",type=float,default=0.0,help="Cutoff by BLAST % identity  (>=),指定値以上を残す")
args = parser.parse_args()



class BLAST2SimMatrix(object):
    def __init__(self, input_file=None, seq_identity_thres=0, edge_weight_thres=0):
        self.blast_file = input_file if input_file else None
        self.seq_identity_thres = seq_identity_thres
        self.edge_weight_thres = edge_weight_thres

    def __repr__(self):
        #blastファイル，パラメータ，とかの情報も吐き出すように。
        if self.entry_set:
            return f"# of Nodes: {len(self.entry_set)}"

    def parse_blast_output(self):
        logger.info("No multi hit BLAST mode! :D")
        logger.info(f"cut off BLAST-hit by % identity < {self.seq_identity_thres}")
        self.blast = defaultdict(dict)
        self.entry_set = set()
        p = re.compile('#')

        in_fh = open(self.blast_file, 'r') if self.blast_file else sys.stdin 
        for line in in_fh:
            if p.match(line):
                continue
            line = line.rstrip()

            query   = line.split('\t')[0]
            subject = line.split('\t')[1]
            score   = float(line.split('\t')[11])
            identity   = float(line.split('\t')[2])

            if identity < self.seq_identity_thres:
                continue
                
            self.blast[query][subject] = score
        
            self.entry_set.add(query)
            self.entry_set.add(subject)
        in_fh.close()
        
    def create_adjacency_matrix(self):
        self.entry_list = list(self.entry_set)
        entry_num = len(self.entry_list)
        self.adj_mat = np.zeros((entry_num, entry_num))
        for i, i_id in enumerate(self.entry_list):
            for j, j_id in enumerate(self.entry_list):
#                ij=0.0; ji=0.0; ii=0.0; jj=0.0;

                try:
                    ij = self.blast[i_id][j_id]
                except KeyError:
                    pass
                try:
                    ji = self.blast[j_id][i_id]
                except KeyError:
                    pass
                try:
                    ii = self.blast[i_id][i_id]
                except KeyError:
                    pass
                try:
                    jj = self.blast[j_id][j_id]
                except KeyError:
                    pass
        
                dist  = ij if ij >= ji else ji
                self_  = ii if ii >= jj else jj 
        
                score = 0.0
                try:
                    score = dist/self_
                except ZeroDivisionError as err: #ないと信じてるが
                    score = dist/1
                    logger.warn("ZeroDivisionError!" + err)
                    othererr.append(err)
        
                self.adj_mat[i][j] = 0.0 if i == j else score #1だと自分自身へのエッジが生まれる

    def output_as_sim_adj_list(self, output_file=None):
        logger.info(f"save edges > {self.edge_weight_thres}")

        out_fh = open(output_file, 'w') if output_file else sys.stdout
        for i, i_id in enumerate(self.entry_list):
            for j, j_id in enumerate(self.entry_list):
                if i >= j and self.adj_mat[i][j] > self.edge_weight_thres:
                    out_fh.write(f"{i_id} {j_id} {self.adj_mat[i][j]}\n") #sim


        logger.info("Output file -> " + output_file) if output_file else logger.info("Output file -> sys.stdout")
        out_fh.close()






if __name__ == "__main__":
    np.set_printoptions(threshold=np.nan)
    start_time = time.time()

    main = BLAST2SimMatrix(input_file=args.input_file, seq_identity_thres=args.cutoff_identity, edge_weight_thres=args.cutoff_edge_weight)
    main.parse_blast_output()
    main.create_adjacency_matrix()
    main.output_as_sim_adj_list()

    logger.debug(main.adj_mat)
    logger.info("Run time: "+ str((time.time() - start_time)/60) + '[m]')
    logger.info("= Run time: "+ str((time.time() - start_time)) + '[s]')
