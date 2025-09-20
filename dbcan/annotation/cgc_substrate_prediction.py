import os
import sys
import logging
import subprocess
from typing import List, Dict
from Bio import SeqIO
import pandas as pd
import numpy as np
import math
import json
import time
from dbcan.configs.cgc_substrate_config import CGCSubstrateConfig
import dbcan.constants.cgc_substrate_prediction_constants as CSP  
# from dbcan.constants import (
#     CAZYME, TC, TF, STP, PUL, NULL, CGC_RESULT_FILE,
#     DBCAN_SUB_OUT_FILE, OVERVIEW_FILE, INPUT_PROTEIN_NAME,
#     CGC_SUB_PREDICTION_FILE, PUL_DIAMOND_FILE, CGC_FAA_FILE,
#     PUL_DIAMOND_DB, PUL_EXCEL_FILE, DBCANPUL_TMP, DBCAN_SUB_TMP,
#     DIAMOND_PUL_EVALUE, CAZYME_FAA_FILE, PUL_FAA_FILE
# )

logger = logging.getLogger(__name__)

def Sum_bitscore(genes):
    return sum([gene.bitscore for gene in genes])

class blastp_hit(object):
    def __init__(self, lines):
        self.qseqid = lines[0]
        self.sseqid = lines[1]
        self.pident = float(lines[2])
        self.length = int(lines[3])
        self.mismatch = int(lines[4])
        self.gapopen  = int(lines[5])
        self.qstart   = int(lines[6])
        self.qend     = int(lines[7])
        self.sstart   = int(lines[8])
        self.send     = int(lines[9])
        self.evalue   = float(lines[10])
        self.bitscore = float(lines[11])
        self.qlen = int(lines[12]) if len(lines) >= 13 else 0
        self.slen = int(lines[13]) if len(lines) >= 14 else 0

    def __repr__(self):
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])

    def format_str(self):
        qseqids = self.qseqid.split("|")
        sseqids = self.sseqid.split(":")
        qtype = qseqids[3] if len(qseqids) > 3 else ""
        stype = sseqids[5] if len(sseqids) > 5 else ""
        if qtype == CSP.CAZYME:
            families = ";".join(qseqids[4:])
            qseqid = qseqids[2] + "|" + qseqids[3] + "|" + families
        else:
            qseqid = qseqids[2] + "|" + qseqids[3]
        if stype == CSP.CAZYME and len(sseqids) > 6:
            sseqid = sseqids[0] + "|" + sseqids[5] + "|" + sseqids[6].replace("|",";")
        else:
            sseqid = sseqids[0] + "|" + (sseqids[5] if len(sseqids) > 5 else "")
        cgcid = qseqids[0] + "|" + qseqids[1]
        pulid = sseqids[1] if len(sseqids) > 1 else ""
        return "\t".join([qseqid,sseqid,cgcid,pulid,str(self.pident),str(self.length),str(self.mismatch),str(self.gapopen),str(self.qstart),str(self.qend),str(self.sstart),str(self.send),str(self.evalue),str(self.bitscore),str(self.qlen),str(self.slen)])

class Gene(object):
    def __init__(self, lines, cluster_number=0):
        self.clusterid = lines[0]
        self.type  =    lines[1]
        self.contig =   lines[2]
        self.Protein_ID = lines[3]
        self.start =    int(lines[4])
        self.end =      int(lines[5])
        self.strand =   lines[6]
        if self.type == CSP.TC:
            self.Protein_Fam = ".".join(lines[7].split("|")[1].split(".")[0:3])
        elif self.type == CSP.NULL:
            self.Protein_Fam = CSP.NULL
        else:
            self.Protein_Fam = lines[7].split("|")[1]
        self.CGC_ID = self.contig + "|" + self.clusterid

    def format_out(self):
        return "\t".join([self.clusterid,self.type,self.contig,self.Protein_ID,str(self.start),str(self.end),self.strand,self.Protein_Fam])

    def __getattr__(self, name):
        return "-"

    def Get_CAzyID(self):
        return self.contig+"|"+self.clusterid+"|"+self.Protein_ID+"|"+self.type+"|"+self.Protein_Fam

class dbCAN_Out(object):
    def __init__(self, filename):
        hits = open(filename).readlines()[1:]
        self.genes: List[Gene] = []
        for line in hits:
            if line.startswith("CGC#"):
                continue
            lines = line.split()
            if not lines:
                continue
            self.genes.append(Gene(lines))

    def __iter__(self):
        return iter(self.genes)

    def ProteinID2genes(self):
        ProteinIDdict = {}
        for gene in self:
            ProteinIDdict[gene.Protein_ID] = gene
        return ProteinIDdict

    def CGCID2genes(self):
        cgcdict = {}
        for gene in self:
            cgcdict.setdefault(gene.CGC_ID, []).append(gene)
        return cgcdict

class CGC(object):
    def __init__(self, genes: List[Gene]):
        self.genes = genes
        self.ID = genes[0].CGC_ID
        self.start = min([gene.start for gene in genes])
        self.end = max([gene.end for gene in genes])
        self.gene_num = len(genes)

    def __iter__(self):
        return iter(self.genes)

class CGC_hub(object):
    def __init__(self, dbcan: dbCAN_Out):
        self.CGCs: List[CGC] = []
        cgcdict = dbcan.CGCID2genes()
        for cgc in cgcdict:
            self.CGCs.append(CGC(cgcdict[cgc]))

    def __iter__(self):
        return iter(self.CGCs)

    def CGCID2CGC(self) -> Dict[str, CGC]:
        return {cgc.ID: cgc for cgc in self}

class dbSub_record(object):
    def __init__(self, lines):
        self.dbcan_sub_subfam = lines[0]
        self.Subfam_Composition = lines[1]
        self.Subfam_EC = lines[2]
        self.Substrate = lines[3] if lines[3] != "-" else ""
        self.hmm_Length = lines[4]
        self.GeneID = lines[5]
        self.GeneLen = lines[6]
        self.E_Value = lines[7]
        self.hmm_Start = lines[8]
        self.hmm_End = lines[9]
        self.Gene_Start = lines[10]
        self.Gene_End = lines[11]
        self.Cov = lines[12]

    def __repr__(self):
        return "\t".join([self.__dict__[name] for name in self.__dict__])

class dbSub(object):
    def __init__(self, filename, dbsub_parameters):
        self.Genes: List[dbSub_record] = []
        if not os.path.exists(filename) or os.path.getsize(filename) == 0:
            return
        with open(filename) as fh:
            next(fh, None)  # skip header
            for line in fh:
                lines = line.rstrip("\n").split("\t")
                if len(lines) < 13:
                    continue
                hmmevalue = float(lines[7])
                hmmcov    = float(lines[12])
                if hmmevalue <= dbsub_parameters.hmmevalue and hmmcov >= dbsub_parameters.hmmcov:
                    self.Genes.append(dbSub_record(lines))

    def __iter__(self):
        return iter(self.Genes)

    def GeneID2gene(self):
        geneid2gene = {}
        for gene in self:
            geneid2gene.setdefault(gene.GeneID, []).append(gene)
        self.geneid2gene = geneid2gene
        return geneid2gene

def replace_black_underline(df):
    col_names = df.columns.tolist()
    for index, value in enumerate(col_names):
        col_names[index] = value.replace(" ", "_")
    df.columns = col_names

class dbCAN_substrate_prediction:
    """
    based on dbCAN-sub output and PUL blast result, predict substrate specificities for each CGC.
    """

    def __init__(self, config: CGCSubstrateConfig):
        self.config = config
        self.input_folder = config.output_dir
        self.cgc_standard_out = os.path.join(self.input_folder, CSP.CGC_RESULT_FILE)
        self.dbsub_out = os.path.join(self.input_folder, CSP.DBCAN_SUB_OUT_FILE)
        self.overview_txt = os.path.join(self.input_folder, CSP.OVERVIEW_FILE)
        self.protein_db = os.path.join(self.input_folder, CSP.INPUT_PROTEIN_NAME)

        self.out = os.path.join(self.input_folder, CSP.CGC_SUB_PREDICTION_FILE)
        self.tmp_folder = self.input_folder
        self.tmp_blastp_out = os.path.join(self.tmp_folder, CSP.PUL_DIAMOND_FILE)
        self.tmp_CAZyme_pep = os.path.join(self.tmp_folder, CSP.CGC_FAA_FILE)

        self.db_dir = config.db_dir
        self.PULdb_diamond = os.path.join(self.db_dir, CSP.PUL_DIAMOND_DB)
        self.pul_excel_filename = os.path.join(self.db_dir, CSP.PUL_EXCEL_FILE)

        self.homologous_parameters  = HitParamter(config)
        self.dbsub_parameters  = dbcan_sub_parameter(config)

        self.odbcan_sub = config.odbcan_sub
        self.odbcanpul = config.odbcanpul
        self.dbcanpul_tmp = CSP.DBCANPUL_TMP
        self.dbcan_sub_tmp = CSP.DBCAN_SUB_TMP

        self.run_dbCAN_sub = True
        self.run_dbCAN_PUL = True

        self.dbcan_sub_CGC2substrates = {}
        self.queryCGC2hit = {}
        self.dbcan_sub_CGC2maxscore = {}

        self.seqs = []
        self.seqid2seq = {}
        self.protid2gene = {}
        self.dbCAN_hits = None
        self.cgcs = None
        self.cgcid2cgc = {}

    def check_input(self):
        ok = True
        if not os.path.exists(self.input_folder):
            logger.error(f"Input folder not found: {self.input_folder}")
            ok = False
        if not os.path.exists(self.cgc_standard_out):
            logger.error(f"CGC standard out not found: {self.cgc_standard_out}")
            ok = False
        if not os.path.exists(self.dbsub_out):
            logger.warning(f"dbsub output not found: {self.dbsub_out}. dbCAN-sub path will be disabled.")
            self.run_dbCAN_sub = False
        if not os.path.exists(self.overview_txt):
            logger.warning(f"Overview file not found: {self.overview_txt}")
        if not os.path.exists(self.protein_db):
            logger.error(f"Protein fasta not found: {self.protein_db}")
            ok = False
        if not os.path.exists(self.PULdb_diamond):
            logger.error(f"PUL DIAMOND DB not found: {self.PULdb_diamond}")
            ok = False
        if not os.path.exists(self.pul_excel_filename):
            logger.error(f"PUL excel not found: {self.pul_excel_filename}")
            ok = False
        return ok

    def extract_seq_in_CGC(self):
        logger.info("Extracting CGC protein sequences...")
        self.dbCAN_hits = dbCAN_Out(self.cgc_standard_out)
        self.seqid2seq = SeqIO.to_dict(SeqIO.parse(self.protein_db, 'fasta'))
        self.protid2gene = self.dbCAN_hits.ProteinID2genes()
        self.seqs = []
        for seqid in self.seqid2seq:
            if seqid in self.protid2gene:
                self.seqid2seq[seqid].id = self.protid2gene[seqid].Get_CAzyID()
                self.seqs.append(self.seqid2seq[seqid])

    def do_blastp_against_dbCANPUL(self):
        logger.info(f"Running DIAMOND against PUL DB: {self.PULdb_diamond}")
        if not self.seqs:
            self.extract_seq_in_CGC()
        os.makedirs(self.tmp_folder, exist_ok=True)
        SeqIO.write(self.seqs, self.tmp_CAZyme_pep, 'fasta')
        threads = self.config.threads if self.config.threads and self.config.threads > 0 else (os.cpu_count() or 1)
        # --outfmt 6 with explicit fields
        outfmt_fields = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen"]
        cmd = [
            "diamond", "blastp", "--max-hsps", "1",
            "--query", self.tmp_CAZyme_pep,
            "--db", self.PULdb_diamond,
            "--evalue", str(CSP.DIAMOND_PUL_EVALUE),
            "--out", self.tmp_blastp_out,
            "--threads", str(threads),
            "--outfmt", "6", *outfmt_fields,
            "--quiet"
        ]
        logger.debug(" ".join(cmd))
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"DIAMOND failed: {e}")
            # write empty output file
            open(self.tmp_blastp_out, 'w').close()

    def read_dbCAN_PUL(self):
        logger.info(f"Reading PUL excel: {self.pul_excel_filename}")
        self.Puls = pd.read_excel(self.pul_excel_filename)
        replace_black_underline(self.Puls)
        if "substrate_final" not in self.Puls.columns:
            raise ValueError("PUL excel missing column 'substrate_final'")

    def read_blastp_result(self, filename):
        logger.info(f"Reading DIAMOND result: {filename}")
        querydict = {}
        if not os.path.exists(filename) or os.path.getsize(filename) == 0:
            self.hitdict = querydict
            return
        with open(filename) as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) < 14:
                    continue
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen = parts[:14]
                qids = qseqid.split("|")
                if len(qids) < 2:
                    logger.warning(f"Invalid qseqid (no '|'): {qseqid}")
                    continue
                queryid = qids[0] + "|" + qids[1]
                if float(evalue) > self.homologous_parameters.evalue_cutoff:
                    continue
                if float(pident) < self.homologous_parameters.identity_cutoff:
                    continue
                if (float(qend) - float(qstart) + 1) / float(qlen) < self.homologous_parameters.coverage_cutoff:
                    continue
                querydict.setdefault(queryid, []).append(blastp_hit(parts))
        self.hitdict = querydict

    def Uniq_blastp_hit(self, blast_list: List[blastp_hit]):
        uniq_sseqid = set()
        genes = []
        homologous_pairs: List[str] = []
        for tmp in blast_list:
            if tmp.sseqid not in uniq_sseqid:
                genes.append(tmp)
                uniq_sseqid.add(tmp.sseqid)
            hit_type = tmp.sseqid.split(":")[-2] if ":" in tmp.sseqid else ""
            qparts = tmp.qseqid.split("|")
            query_type = qparts[3] if len(qparts) > 3 else ""
            if hit_type and query_type:
                homologous_pairs.append(query_type + "-" + hit_type)
        CAZyme_pairs_num = homologous_pairs.count("CAZyme-CAZyme")
        cond_base = (
            len(uniq_sseqid) >= self.homologous_parameters.uqcgn and
            len({h.qseqid for h in blast_list}) >= self.homologous_parameters.upghn and
            CAZyme_pairs_num >= self.homologous_parameters.cpn and
            len(homologous_pairs) >= self.homologous_parameters.tpn
        )
        if not cond_base:
            return -1, homologous_pairs
        if not self.homologous_parameters.ept:
            score = Sum_bitscore(genes)
            return (score, homologous_pairs) if (score / max(1, len(uniq_sseqid)) >= self.homologous_parameters.bitscore_cutoff) else (-1, homologous_pairs)
        # extra signature pairs needed
        signature_pairs = self.homologous_parameters.ept
        signature_pairs_num = self.homologous_parameters.eptn
        pair_signal = 0
        for i, signature_pair in enumerate(signature_pairs):
            if homologous_pairs.count(signature_pair) >= int(signature_pairs_num[i]):
                pair_signal += 1
        if pair_signal == len(signature_pairs):
            score = Sum_bitscore(genes)
            return (score, homologous_pairs) if (score / max(1, len(uniq_sseqid)) >= self.homologous_parameters.bitscore_cutoff) else (-1, homologous_pairs)
        return -1, homologous_pairs

    def dbcan_sub_read_cgc(self):
        if not self.cgcid2cgc:
            self.dbCAN_hits = dbCAN_Out(self.cgc_standard_out)
            self.cgcs = CGC_hub(self.dbCAN_hits)
            self.cgcid2cgc = self.cgcs.CGCID2CGC()

    def analyze_blastp_out(self):
        logger.info("Analyzing DIAMOND result...")
        if not os.path.exists(self.tmp_blastp_out) or os.path.getsize(self.tmp_blastp_out) == 0:
            logger.warning(f"BLASTP output empty: {self.tmp_blastp_out}")
            os.makedirs(os.path.join(self.input_folder, "synteny_pdf"), exist_ok=True)
            with open(self.out, 'w') as f:
                f.write("#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignature pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score\n")
            # write empty output file
            self.hitdict = {}
            self.queryCGC2scores = {}
            self.queryCGC2pulids = {}
            self.queryCGCmapedtypes = {}
            self.queryCGC2blastp_hits = {}
            return
        self.read_blastp_result(self.tmp_blastp_out)
        self.cgcs = CGC_hub(self.dbCAN_hits)
        self.cgcid2cgc = self.cgcs.CGCID2CGC()

        queryCGC2scores = {}; queryCGC2pulids = {}; queryCGC2Mapped_types = {}
        self.queryCGC2blastp_hits = {}
        for hit in self.hitdict:
            tmp_dict = {}
            for hit2 in self.hitdict[hit]:
                pulid = hit2.sseqid.split("_")[0] if "_" in hit2.sseqid else hit2.sseqid
                tmp_dict.setdefault(pulid, []).append(hit2)
            scores = []; pulids = []; maped_types = []
            for pulid in tmp_dict:
                score, homologous_pairs = self.Uniq_blastp_hit(tmp_dict[pulid])
                if score > 0:
                    scores.append(score)
                    pulids.append(pulid)
                    maped_types.append(homologous_pairs)
                    self.queryCGC2blastp_hits.setdefault(hit, []).append(tmp_dict[pulid])
            queryCGC2scores[hit] = scores
            queryCGC2pulids[hit] = pulids
            queryCGC2Mapped_types[hit] = maped_types
        self.queryCGC2scores = queryCGC2scores
        self.queryCGC2pulids = queryCGC2pulids
        self.queryCGCmapedtypes = queryCGC2Mapped_types

    def get_best_pul_hit_and_blastphit(self):
        logger.info("Selecting best PUL hit per CGC...")
        self.read_dbCAN_PUL()
        queryCGC2hit = {}
        self.queryCGC_best_genes_blastp_hit = {}
        self.queryCGC_CGChits_genes_blastp_hit = {}
        for queryCGC in self.queryCGC2scores:
            scores = self.queryCGC2scores[queryCGC]
            if len(scores) == 0:
                continue
            score_orders = np.argsort(scores)
            max_score_index = score_orders[-1]
            best_score = scores[max_score_index]
            bestpulid = self.queryCGC2pulids[queryCGC][max_score_index]
            substrate = self.Puls[self.Puls["ID"] == bestpulid].substrate_final.values[0]
            mapped_types = self.queryCGCmapedtypes[queryCGC][max_score_index]
            queryCGC2hit[queryCGC] = PULhit(best_score, bestpulid, substrate, mapped_types)
            self.queryCGC_best_genes_blastp_hit[queryCGC] = self.queryCGC2blastp_hits.get(queryCGC, [])[max_score_index] if self.queryCGC2blastp_hits.get(queryCGC) else []
            for idx in score_orders:
                if self.queryCGC2blastp_hits.get(queryCGC):
                    self.queryCGC_CGChits_genes_blastp_hit.setdefault(queryCGC, []).append(self.queryCGC2blastp_hits[queryCGC][idx])
        self.queryCGC2hit = queryCGC2hit

    def Read_CAZyme_substrate(self):
        logger.info(f"Reading dbsub outfile: {self.dbsub_out}")
        self.CAZyme2substrate = dbSub(self.dbsub_out, self.dbsub_parameters)
        self.geneid2dbsub = self.CAZyme2substrate.GeneID2gene() if getattr(self.CAZyme2substrate, "Genes", None) else {}
        cgcid2sub = {}
        for cgcid in self.cgcid2cgc:
            for gene in self.cgcid2cgc[cgcid]:
                if gene.Protein_ID in self.geneid2dbsub:
                    cgcid2sub.setdefault(cgcid, []).append(self.geneid2dbsub[gene.Protein_ID])

    def dbcan_sub_subfamily_substrate_prediction(self):
        self.Read_CAZyme_substrate()
        self.dbcan_sub_read_cgc()
        self.CGC2substrate_dbcan_sub()
        self.substrate_scoring_dbcan_sub()

    def substrate_scoring_dbcan_sub(self):
        logger.info("Scoring dbCAN-sub substrates...")
        finalsub = {}; finalscores = {}; finalranks = {}; finalmaxscore = {}
        cgcid2sub = self.cgcid2substrate_dbcan_sub
        for cgcid in cgcid2sub:
            if self.cgcid2CAZyme_substrate_num[cgcid] < self.dbsub_parameters.num_of_domains_substrate_cutoff:
                continue
            if self.cgcid2CAZyme_domain_substrate_num[cgcid] < self.dbsub_parameters.num_of_protein_shared_substrate_cutoff:
                continue
            scores = {}; ranks = []
            for subs_list in cgcid2sub[cgcid]:
                for subs in subs_list:
                    subs_str = subs.Substrate.replace("and", ",").replace(" ", "")
                    for tmp_sub in set(filter(None, subs_str.split(","))):
                        for i, alt in enumerate(tmp_sub.split("|")):
                            scores.setdefault(alt, []).append(math.pow(2, -i))
            for sub in scores:
                scores[sub] = sum(scores[sub])
                ranks.append(f"{sub}:{scores[sub]}")
            if not scores:
                continue
            finalscores[cgcid] = scores
            max_score = max(scores.values())
            if max_score < self.dbsub_parameters.dbcan_substrate_scors:
                continue
            finalmaxscore[cgcid] = max_score
            final_subs = [sub for sub, sc in scores.items() if sc == max_score]
            finalsub[cgcid] = ",".join(final_subs)
            finalranks[cgcid] = ranks

        self.dbcan_sub_CGC2substrates = finalsub
        self.dbcan_sub_CGC2scores = finalscores
        self.dbcan_sub_substrate_score = finalranks
        self.dbcan_sub_CGC2maxscore = finalmaxscore

    def CGC2substrate_dbcan_sub(self):
        cgcid2sub = {}
        for cgcid in self.cgcid2cgc:
            for gene in self.cgcid2cgc[cgcid]:
                if getattr(self, "geneid2dbsub", None) and gene.Protein_ID in self.geneid2dbsub:
                    cgcid2sub.setdefault(cgcid, []).append(self.geneid2dbsub[gene.Protein_ID])
        self.cgcid2substrate_dbcan_sub = cgcid2sub

        self.cgcid2CAZyme_domain_substrate_num = {}
        for cgcid in cgcid2sub:
            gene_ids = []
            for dbsub_records in cgcid2sub[cgcid]:
                for dbsub_record in dbsub_records:
                    if dbsub_record.Substrate:
                        gene_ids.append(dbsub_record.GeneID)
            self.cgcid2CAZyme_domain_substrate_num[cgcid] = len(set(gene_ids))

        cgcid2substrate_CAZyme_num = {}
        for cgcid in cgcid2sub:
            cnt = 0
            for dbsub_records in cgcid2sub[cgcid]:
                for dbsub_record in dbsub_records:
                    if dbsub_record.Substrate:
                        cnt += 1
            cgcid2substrate_CAZyme_num[cgcid] = cnt
        self.cgcid2CAZyme_substrate_num = cgcid2substrate_CAZyme_num

    def substrate_predict(self):
        if self.run_dbCAN_PUL:
            self.dbCAN_PUL_substrate_predict()
        if self.run_dbCAN_sub:
            self.dbcan_sub_subfamily_substrate_prediction()

    def dbCAN_PUL_substrate_predict(self):
        self.extract_seq_in_CGC()
        self.do_blastp_against_dbCANPUL()
        self.analyze_blastp_out()
        self.get_best_pul_hit_and_blastphit()

    def result_print_to_file(self):
        shared_cgcids = self.queryCGC2hit.keys() | self.dbcan_sub_CGC2substrates.keys()
        logger.info(f"Writing substrate prediction result to file: {self.out}")
        with open(self.out, 'w') as f:
            f.write("#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignature pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score\n")
            for cgcid in shared_cgcids:
                dbcan_pul_part = self.queryCGC2hit.get(cgcid, "")
                dbcan_sub_substate = self.dbcan_sub_CGC2substrates.get(cgcid, "")
                PULID = dbcan_pul_part.pulid if dbcan_pul_part else ""
                dbcan_pul_sub = dbcan_pul_part.substrate if dbcan_pul_part else ""
                bitscore = dbcan_pul_part.score if dbcan_pul_part else ""
                sig_pairs = ";".join(dbcan_pul_part.maped_types) if dbcan_pul_part else ""
                dbcan_sub_maxscore = self.dbcan_sub_CGC2maxscore.get(cgcid, "")
                f.write(f"{cgcid}\t{PULID}\t{dbcan_pul_sub}\t{bitscore}\t{sig_pairs}\t{dbcan_sub_substate}\t{dbcan_sub_maxscore}\n")

class PULhit(object):
    def __init__(self, score, pulid, substrate, mapped_types):
        self.score = score
        self.substrate = substrate
        self.maped_types = mapped_types
        self.pulid = pulid
    def __repr__(self):
        return "\t".join([self.pulid, str(self.score), self.substrate, ",".join(self.maped_types)])

def cgc_substrate_prediction(config: CGCSubstrateConfig):
    sp = dbCAN_substrate_prediction(config)
    if not sp.check_input():
        # write empty result and return
        os.makedirs(config.output_dir, exist_ok=True)
        out = os.path.join(config.output_dir, CSP.CGC_SUB_PREDICTION_FILE)
        with open(out, 'w') as f:
            f.write("#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignature pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score\n")
        logger.error("Required input missing. Wrote empty substrate prediction file.")
        return
    t0 = time.time()
    sp.substrate_predict()
    t1 = time.time()
    logger.info(f"Substrate prediction done! {t1 - t0:.2f}s")
    sp.result_print_to_file()
    if sp.odbcan_sub:
        sp.dbcan_sub_intermediate_file()
class HitParamter(object):
    '''
    design for parameters, how to identify real homologous hit for genes in GCG
    '''
    def __init__(self,config):
        self.config = config

        self.upghn = config.uniq_pul_gene_hit_num
        self.uqcgn = config.uniq_query_cgc_gene_num
        self.cpn = config.CAZyme_pair_num
        self.tpn = config.total_pair_num
        self.ept = config.extra_pair_type.split(",") if config.extra_pair_type else None
        self.eptn = config.extra_pair_type_num.split(",") if self.ept else 0
        self.identity_cutoff = config.identity_cutoff
        self.coverage_cutoff  = config.coverage_cutoff
        self.bitscore_cutoff = config.bitscore_cutoff
        self.evalue_cutoff = config.evalue_cutoff

        ### check the additional requires for the other signature pairs
        if self.ept and len(self.ept) != len(self.eptn):
            print(f"The optional chocices of {self.ept} is not equal to {self.eptn}.",file=sys.stderr)
            exit()

    def __repr__(self):
        return "\n".join([name + ": " +str(self.__dict__[name]) for name in self.__dict__])

class dbcan_sub_parameter(object):
    def __init__(self,config):
        self.config = config
        self.hmmevalue = config.hmmevalue
        self.hmmcov    = config.hmmcov
        self.num_of_protein_shared_substrate_cutoff = config.num_of_protein_substrate_cutoff
        self.num_of_domains_substrate_cutoff = config.num_of_domains_substrate_cutoff
        self.dbcan_substrate_scors =  config.substrate_scors

    def __repr__(self):
        return "\n".join([name + ": " +str(self.__dict__[name]) for name in self.__dict__])
def clean_sub(sub):
    subs = sub.Substrate
    subs = subs.replace("and",",") ### based some naming habit by someone
    subss = subs.split(",")
    subss = set(subss)
    tmp_subs = []
    for tmp_sub in subss: ### loop for substrate
        if not tmp_sub: ### exclude "" substrate come from "-"
            continue
        tmp_sub = tmp_sub.strip(" ") ### remove blank in the two ends of substrate,
        tmp_subs.extend(tmp_sub.split("|")) ### some substrates combined by "|"
    return list(set(tmp_subs))

def clean_EC(sub):
    subs = sub.Subfam_EC
    tmp_subs = []
    for tmp_sub in subs.split("|"): ### loop for substrate
        if not tmp_sub or tmp_sub =="-": ### exclude "" substrate come from "-"
            continue
        tmp_sub = tmp_sub.strip(" ") ### remove blank in the two ends of substrate,
        tmp_subs.append(tmp_sub.split(":")[0])
    return list(set(tmp_subs))