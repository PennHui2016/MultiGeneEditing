import sys
import scipy.io as sio
import csv
import numpy as np
import math
import pandas as pd
import sys
import multiprocessing
import os
import copy
import itertools as iter
from functools import reduce
import statistics


#rootdir='../'
rootdir=os.path.abspath(os.path.join(os.path.dirname(__file__),'../')).replace('\\','/')+'/'

def obtain_gene_info():
    file=rootdir+'data/cancer_gene_info1.csv'
    ids=[]
    symbols=[]
    chrs=[]
    strands=[]
    starts=[]
    ends=[]
    f=open(file)
    for line in f:
        s=line.split('\t')
        ids.append(s[0])
        symbols.append(s[1])
        chrs.append(s[2])
        starts.append(int(s[3]))
        ends.append(int(s[4]))
        strands.append(s[5])

    return ids, symbols, chrs, starts, ends, strands

def readOnTars(gene_name):
    dir_ontar = rootdir + 'data/predict_genes_results/'
    file_name=dir_ontar+gene_name+'.csv'
    f=open(file_name,'r')
    fdict_reader = csv.DictReader(f)
    sgRNAs_ex=[]
    spacers=[]
    strand=[]
    poss=[]
    pre_on_scores=[]
    for line in fdict_reader:
        sgRNAs_ex.append(line['spacer_ex'])
        spacer=line['spacer']
        spacers.append(spacer)
        strand.append(line['strand'])
        poss.append(line['position'])
        pre_on_scores.append(line['TSAM_predict'])
    return  sgRNAs_ex, spacers, strand, poss,pre_on_scores

def reverse_seq(seq):
    re_seq=''
    for i in range(0,len(seq)):
        if seq[i]=='A':
            re_seq=re_seq+'T'
        elif seq[i]=='G':
            re_seq=re_seq+'C'
        elif seq[i]=='C':
            re_seq=re_seq+'G'
        elif seq[i]=='T':
            re_seq=re_seq+'A'

    return re_seq

def mutate_score(mutate_seq, ref_seq, mutation_types, mutation_scores,mutation_position, ref_score):
    ref_seq=ref_seq.replace('T','U')
    mutate_seq=reverse_seq(mutate_seq)
    mu_score=float(ref_score)
    for i in range(0,len(mutate_seq)):
        index=[j for j, X in enumerate(mutation_position) if X == str(i+1)]
        mutates=list(np.array(mutation_types)[np.array(index)])
        mutates_scores=list(np.array(mutation_scores)[np.array(index)])
        mutate='r'+ref_seq[i]+':'+'d'+mutate_seq[i]
        if mutate in mutates:
            ind=mutates.index(mutate)
            mu_score=mu_score*float(mutates_scores[ind])
    return mu_score

def obtain_exons(file):
    gene_chr=[]
    gene_start=[]
    gene_end=[]
    gene_strand=[]
    gene_id=[]
    gene_name=[]

    with open(file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                lin = line.strip().split('\t')
                chr=lin[0]
                Type = lin[2]
                if Type == "gene":  ##extract exon information
                    info = lin[8].split('; ')
                    gene_chr.append(chr)
                    gene_start.append(int(lin[3]))
                    gene_end.append(int(lin[4]))
                    gene_strand.append(lin[6])
                    for str in info:
                        sstr = str.split(" ")
                        if sstr[0] == "gene_id":
                            gene_id.append(sstr[1].replace('"', ''))
                        elif sstr[0] == "gene_name":
                            gene_name.append(sstr[1].replace('"', ''))

    return gene_chr,gene_start,gene_end,gene_strand,gene_id,gene_name

mutation_types=[]
mutation_scores=[]
mutation_position=[]
f=open(rootdir+'data/mismatch_efficiency.txt')
for line in f:
    contents=line.split('\t')
    mutation_types.append(contents[0])
    mutation_position.append(contents[1])
    mutation_scores.append(float(contents[5]))

# def TGE(gene1, gene2, gene_chr,gene_start,gene_end,gene_strand,gene_id,gene_name):
#     off_dir = '/projects/EGANBT/pre_cancer_offs/'
#     gene_all = pd.DataFrame(zip(gene_id, gene_name, gene_chr, gene_start, gene_end), columns=['id', 'name', 'chr', 'start', 'end'])
#     index_gene1 = np.where(gene_all['id'] == gene1)
#     index_gene2 = np.where(gene_all['id'] == gene2)
#     chr_g1 = gene_all['chr'][index_gene1[0][0]]
#     chr_g2 = gene_all['chr'][index_gene2[0][0]]
#     gene1_start=gene_all['start'][index_gene1[0][0]]
#     gene1_end = gene_all['end'][index_gene1[0][0]]
#     gene2_start = gene_all['start'][index_gene2[0][0]]
#     gene2_end = gene_all['end'][index_gene2[0][0]]
#     dir = '/projects/EGANBT/multi_gene_editing/single_editing_modify/'
#     file_path1 = dir + gene1 + '_' + str(0.5) + '.csv'
#     sgs1 = pd.read_csv(file_path1)
#     file_path2 = dir + gene2 + '_' + str(0.5) + '.csv'
#     sgs2 = pd.read_csv(file_path2)
#
#     spacer1=[]
#     spacer_ex1=[]
#     chr1=[]
#     pos1=[]
#     strand1=[]
#     genes1=[]
#
#     spacer2 = []
#     spacer_ex2 = []
#     chr2 = []
#     pos2 = []
#     strand2 = []
#     genes2=[]
#
#     select_type = []
#
#     ce1=[]
#     ce2=[]
#     ce_avg=[]
#
#     off_num1=[]
#     off_num2=[]
#     total_off=[]
#
#     diff1=[]
#     diff2=[]
#     diff_avg=[]
#
#
#
#     common_sp=np.intersect1d(sgs1['spacer'],sgs2['spacer'])
#     if len(common_sp)>0:
#         for k in common_sp:
#             ind1=sgs1['spacer'].tolist().index(k)
#             ind2=sgs2['spacer'].tolist().index(k)
#             cut1=float(sgs1['on_score'][ind1][1:len(sgs1['on_score'][ind1])-1])
#             cut2 = float(sgs2['on_score'][ind2][1:len(sgs2['on_score'][ind2]) - 1])
#             cut_avg=2*cut1*cut2/(cut1+cut2)
#
#             off_num=sgs1['num_off'][ind1]
#
#             on_diff1=float(sgs1['mean_cut_diff'][ind1])
#             on_diff2 = float(sgs2['mean_cut_diff'][ind2])
#
#             diff_avg0=2*on_diff1*on_diff2/(on_diff1+on_diff2)
#
#             spacer1.append(k)
#             spacer_ex1.append(sgs1['spacer_ex'][ind1])
#             chr1.append(chr_g1)
#             pos1.append(sgs1['position'][ind1])
#             strand1.append(sgs1['strand'][ind1])
#             genes1.append(gene1)
#
#             spacer2.append(k)
#             spacer_ex2.append(sgs2['spacer_ex'][ind2])
#             chr2.append(chr_g2)
#             pos2.append(sgs2['position'][ind2])
#             strand2.append(sgs2['strand'][ind2])
#             genes2.append(gene2)
#
#             select_type.append('common spacer')
#
#             ce1.append(cut1)
#             ce2.append(cut2)
#             ce_avg.append(cut_avg)
#
#             off_num1.append(sgs1['num_off'][ind1])
#             off_num2.append(sgs1['num_off'][ind1])
#             total_off.append(sgs1['num_off'][ind1])
#
#             diff1.append(on_diff1)
#             diff2.append(on_diff2)
#             diff_avg.append(diff_avg0)
#
#
#     for i in range(0, math.floor(len(sgs1)*0.2)):
#         strand_g1 = sgs1['strand'][i]
#         pos_g1 = sgs1['position'][i]
#         cut_diff_g1 = sgs1['mean_cut_diff'][i]
#         ref_seq = sgs1['spacer'][i]
#         ref_seq_ex = sgs1['spacer_ex'][i]
#         on_score = float(sgs1['on_score'][i][1:len(sgs1['on_score'][i]) - 1])
#         off_file1 = off_dir + gene1 + '_' + sgs1['spacer_ex'][i] + '.csv'
#         offs1 = pd.read_csv(off_file1, header=None, sep=' ', names=['spacer', 'chr', 'position', 'strand', 'misnum', 'off_score'])
#         dfs1 = offs1[(offs1['chr'] == chr_g2) & (offs1['position'] >= gene2_start) & (offs1['position'] <= gene2_end)]
#         if len(dfs1)>0:
#             off_cut_scores = []
#
#             for j in range(0, len(offs1)):
#                 mut_seq = offs1['spacer'][j][0:20]
#
#                 mu_score = mutate_score(mut_seq, ref_seq, mutation_types, mutation_scores, mutation_position, on_score)
#                 off_cut_scores.append(mu_score)
#
#             if len(off_cut_scores) <= 1:
#                 mean_diff = [sum(off_cut_scores)]
#             else:
#                 mean_diff = [(i * (len(off_cut_scores) - 1) - (sum(off_cut_scores) - i)) / (len(off_cut_scores) - 1) * (
#                             (sum(j < i for j in off_cut_scores) - 1) / len(off_cut_scores)) for i in off_cut_scores]
#
#             off_sps=dfs1['spacer'].tolist()
#             ids=[offs1['spacer'].tolist().index(t) for t in off_sps]
#             ids=np.asarray(ids)
#
#             spacer1.extend([ref_seq for t in ids])
#             spacer_ex1.extend([ref_seq_ex for t in ids])
#             chr1.extend([chr_g1 for t in ids])
#             pos1.extend([pos_g1 for t in ids])
#             strand1.extend([strand_g1 for t in ids])
#             genes1.extend([gene1 for t in ids])
#
#             spacer2.extend(offs1['spacer'][ids])
#             spacer_ex2.extend(offs1['spacer'][ids])
#             chr2.extend([chr_g2 for t in ids])
#             pos2.extend(offs1['position'][ids])
#             strand2.extend(offs1['strand'][ids])
#             genes2.extend([gene2 for t in ids])
#
#             select_type.extend(['off target' for t in ids])
#
#             ce1.extend([on_score for t in ids])
#             ce2.extend(list(np.array(off_cut_scores)[ids]))
#             ce_avg.extend([2*on_score*off_cut_scores[t]/(on_score+off_cut_scores[t]) for t in ids])
#
#             off_num1.extend([sgs1['num_off'][i]-1 for t in ids])
#             off_num2.extend([sgs1['num_off'][i]-1 for t in ids])
#             total_off.extend([sgs1['num_off'][i]-1 for t in ids])
#
#             diff1.extend([cut_diff_g1 for t in ids])
#             diff2.extend(list(np.array(mean_diff)[ids]))
#             diff_avg.extend([2*cut_diff_g1*mean_diff[t]/(cut_diff_g1+mean_diff[t]) for t in ids])
#
#
#
#
#     for j in range(0, math.floor(len(sgs2) * 0.2)):
#         strand_g2 = sgs2['strand'][j]
#         pos_g2 = sgs2['position'][j]
#         cut_diff_g2 = sgs2['mean_cut_diff'][j]
#         ref_seq = sgs2['spacer'][j]
#         ref_seq_ex = sgs2['spacer_ex'][j]
#         on_score = float(sgs2['on_score'][j][1:len(sgs2['on_score'][j]) - 1])
#         off_file2 = off_dir + gene2 + '_' + sgs2['spacer_ex'][j] + '.csv'
#         offs2 = pd.read_csv(off_file2, header=None, sep=' ',
#                             names=['spacer', 'chr', 'position', 'strand', 'misnum', 'off_score'])
#         dfs2 = offs2[(offs2['chr'] == chr_g1) & (offs2['position'] >= gene1_start) & (offs2['position'] <= gene1_end)]
#         if len(dfs2) > 0:
#             off_cut_scores = []
#
#             for n in range(0, len(offs2)):
#                 mut_seq = offs2['spacer'][n][0:20]
#
#                 mu_score = mutate_score(mut_seq, ref_seq, mutation_types, mutation_scores, mutation_position, on_score)
#                 off_cut_scores.append(mu_score)
#
#             if len(off_cut_scores) <= 1:
#                 mean_diff = [sum(off_cut_scores)]
#             else:
#                 mean_diff = [(i * (len(off_cut_scores) - 1) - (sum(off_cut_scores) - i)) / (len(off_cut_scores) - 1) * (
#                         (sum(j < i for j in off_cut_scores) - 1) / len(off_cut_scores)) for i in off_cut_scores]
#
#             off_sps = dfs2['spacer'].tolist()
#             ids = [offs2['spacer'].tolist().index(t) for t in off_sps]
#             ids = np.asarray(ids)
#
#             spacer1.extend([ref_seq for t in ids])
#             spacer_ex1.extend([ref_seq_ex for t in ids])
#             chr1.extend([chr_g2 for t in ids])
#             pos1.extend([pos_g2 for t in ids])
#             strand1.extend([strand_g2 for t in ids])
#             genes1.extend([gene2 for t in ids])
#
#             spacer2.extend(offs2['spacer'][ids])
#             spacer_ex2.extend(offs2['spacer'][ids])
#             chr2.extend([chr_g1 for t in ids])
#             pos2.extend(offs2['position'][ids])
#             strand2.extend(offs2['strand'][ids])
#             genes2.extend([gene1 for t in ids])
#
#             select_type.extend(['off target' for t in ids])
#
#             ce1.extend([on_score for t in ids])
#             ce2.extend(list(np.array(off_cut_scores)[ids]))
#             ce_avg.extend([2 * on_score * off_cut_scores[t] / (on_score + off_cut_scores[t]) for t in ids])
#
#             off_num1.extend([sgs2['num_off'][j] - 1 for t in ids])
#             off_num2.extend([sgs2['num_off'][j] - 1 for t in ids])
#             total_off.extend([sgs2['num_off'][j] - 1 for t in ids])
#
#             diff1.extend([cut_diff_g2 for t in ids])
#             diff2.extend(list(np.array(mean_diff)[ids]))
#             diff_avg.extend([2 * cut_diff_g2 * mean_diff[t] / (cut_diff_g2 + mean_diff[t]) for t in ids])
#
#     spacer1.append(sgs1['spacer'][0])
#     spacer_ex1.append(sgs1['spacer_ex'][0])
#     chr1.append(chr_g1)
#     pos1.append(sgs1['position'][0])
#     strand1.append(sgs1['strand'][0])
#     genes1.append(gene1)
#
#     spacer2.append(sgs2['spacer'][0])
#     spacer_ex2.append(sgs2['spacer_ex'][0])
#     chr2.append(chr_g2)
#     pos2.append(sgs2['position'][0])
#     strand2.append(sgs2['strand'][0])
#     genes2.append(gene2)
#
#     select_type.append('independent design')
#
#     ce1.append(float(sgs1['on_score'][0][1:len(sgs1['on_score'][0]) - 1]))
#     ce2.append(float(sgs2['on_score'][0][1:len(sgs2['on_score'][0]) - 1]))
#     ce_avg.append(2 * float(sgs1['on_score'][0][1:len(sgs1['on_score'][0]) - 1]) * float(sgs2['on_score'][0][1:len(sgs2['on_score'][0]) - 1]) / (float(sgs1['on_score'][0][1:len(sgs1['on_score'][0]) - 1]) + float(sgs2['on_score'][0][1:len(sgs2['on_score'][0]) - 1])))
#
#     off_num1.append(sgs1['num_off'][0])
#     off_num2.append(sgs2['num_off'][0])
#     total_off.append(sgs1['num_off'][0]+sgs2['num_off'][0])
#
#     diff1.append(sgs1['mean_cut_diff'][0])
#     diff2.append(sgs2['mean_cut_diff'][0])
#     diff_avg.append(2 * sgs1['mean_cut_diff'][0] * sgs2['mean_cut_diff'][0] / (sgs1['mean_cut_diff'][0]+sgs2['mean_cut_diff'][0]))
#
#     rank_avg_cut = ranks(ce_avg, True)
#     rank_off = ranks(total_off, False)
#     rank_diff = ranks(diff_avg, True)
#     final_rank = []
#     for s in range(0, len(ce_avg)):
#         # rank = rank_avg_cut[k][2] * ratio + (1 - ratio) * rank_diff[k][2]
#         # rank = (rank_avg_cut[k][2] * 0.5 + 0.5 * rank_off[k][2])*ratio+rank_diff[k][2]*(1-ratio)
#         rank = (rank_avg_cut[s][2] + rank_off[s][2] + rank_diff[s][2]) / 3
#         final_rank.append(rank)
#
#     final_index = np.zeros((len(final_rank), 2))
#     for i in range(0, len(final_index)):
#         final_index[i, 0] = final_rank[i]
#         final_index[i, 1] = i
#     # Out = np.hstack((final_rank, final_index))
#     out = sorted(final_index, key=lambda row: row[0], reverse=False)
#     file_path = '/projects/EGANBT/multi_gene_editing/two_gene_editing/' + gene1 + '_' + gene2 + '.csv'
#     with open(file_path, 'w') as csv_file:
#         fieldnames = ['select_type','spacer1', 'spacer1_ex', 'spacer2', 'spacer2_ex', 'gene1', 'gene1_chr', 'gene1_strand',
#                       'gene1_position', 'gene2',
#                       'gene2_chr', 'gene2_strand', 'gene2_position', 'cut_eff_g1', 'cut_eff_g2', 'average_cut_on',
#                       'rank_on', 'off_num1', 'off_num2',
#                       'total_off-tar_num', 'rank_off', 'diff1', 'diff2', 'mean_diff', 'rank_diff']
#         writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
#         writer.writeheader()
#         for i in range(0, len(ce_avg)):
#             ind = int(out[i][1])
#             writer.writerow(
#                 {'select_type':select_type[ind],'spacer1': spacer1[ind], 'spacer1_ex': spacer_ex1[ind], 'spacer2': spacer2[ind], 'spacer2_ex':spacer_ex2[ind],
#                  'gene1': genes1[ind], 'gene1_chr': chr1[ind],
#                  'gene1_strand': strand1[ind], 'gene1_position': pos1[ind],
#                  'gene2': genes2[ind], 'gene2_chr': chr2[ind], 'gene2_strand': strand2[ind],
#                  'gene2_position': pos2[ind],
#                  'cut_eff_g1': ce1[ind], 'cut_eff_g2': ce2[ind],
#                  'average_cut_on': ce_avg[ind], 'rank_on': rank_avg_cut[ind][2],'off_num1':off_num1[ind],'off_num2':off_num2[ind],
#                  'total_off-tar_num': total_off[ind], "rank_off": rank_off[ind][2],
#                  'diff1': diff1[ind], 'diff2': diff2[ind],
#                  'mean_diff': diff_avg[ind],
#                  'rank_diff': rank_diff[ind][2]})
#     return file_path

def ranks(score,Type):
    score_index = np.zeros((len(score), 2))
    index=np.zeros((len(score), 1))
    for i in range(0, len(score)):
        score_index[i, 0] = score[i]
        score_index[i,1] = i
        index[i,0] = i
    #Out = np.hstack((score, score_index))
    out = sorted(score_index, key=lambda row: row[0], reverse=Type)
    out_index=np.hstack((out, index))

    out_final = sorted(out_index, key=lambda row: row[1], reverse=False)
    return out_final

# def ranks_combined(select_spacer, candi_pair_1, candi_pair_2, candi_onTs, candi_pair_offs_nums, candi_type, candi_pair_pos1, candi_pair_pos2, candi_pair_stra1, candi_pair_stra2, ratios):
#     on_scores=candi_onTs[2:len(candi_onTs[:,0]),2]
#     offT_nums=candi_pair_offs_nums[2:len(candi_pair_offs_nums[:,0]),2]
#     rank_on=ranks(on_scores,True)
#     rank_off=ranks(offT_nums,False)
#     final_rank = []
#     for i in range(0, len(on_scores)):
#         rank = rank_on[i][2] * ratios + (1 - ratios) * rank_off[i][2]
#         final_rank.append(rank)
#
#     final_index = np.zeros((len(final_rank), 2))
#     for i in range(0, len(final_index)):
#         final_index[i, 0] = final_rank[i]
#         final_index[i, 1] = i
#     # Out = np.hstack((final_rank, final_index))
#     out = sorted(final_index, key=lambda row: row[0], reverse=False)
#     file_path = rootdir+'predict_results/two_gene_overall_ranks_1_6.csv'
#     with open(file_path, 'w') as csv_file:
#         fieldnames = ['select_spacer', 'select_type','average_cut_on', 'rank_on', 'total_off-tar_num', 'rank_off', 'on_tar_seq_1', 'position_1', "strand_1", "on_score_1", "num_off_1",
#                       'on_tar_seq_2', 'position_2', "strand_2", "on_score_2", "num_off_2"]
#         writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
#         writer.writeheader()
#         for i in range(0, len(on_scores)):
#             ind = int(out[i][1])
#             writer.writerow(
#                 {'select_spacer': select_spacer[ind], 'select_type':candi_type[ind], 'average_cut_on': candi_onTs[ind+2,2], 'rank_on': rank_on[ind][2],
#                   'total_off-tar_num':candi_pair_offs_nums[ind+2,2], "rank_off": rank_off[ind][2], 'on_tar_seq_1':candi_pair_1[ind],
#                  'position_1':candi_pair_pos1[ind], "strand_1":candi_pair_stra1[ind], "on_score_1":candi_onTs[ind+2,0], "num_off_1":candi_pair_offs_nums[ind+2,0],
#                  'on_tar_seq_2':candi_pair_2[ind], 'position_2':candi_pair_pos2[ind], "strand_2":candi_pair_stra2[ind],
#                  "on_score_2":candi_onTs[ind+2,1], "num_off_2":candi_pair_offs_nums[ind+2,1]})
#     return file_path

def gene_location(id,gene_all):
    #off_dir = '/projects/EGANBT/pre_cancer_offs/'

    index_gene = np.where(gene_all['id'] == id)
    chr_g = gene_all['chr'][index_gene[0][0]]
    gene_start=gene_all['start'][index_gene[0][0]]
    gene_end = gene_all['end'][index_gene[0][0]]
    dir = rootdir+'data/single_editing_modify/'
    file_path = dir + id + '_' + str(0.5) + '.csv'
    sgs = pd.read_csv(file_path)
    return chr_g, gene_start, gene_end, sgs

def comine_pairs(sgss, dict_pair,selects,genes,chr_gs):
    types=[]
    sgRNAs=[]
    gene_alls=[]
    chrs=[]
    strands=[]
    poss=[]
    on_effs=[]
    off_nums=[]
    cut_diffs=[]
    if len(genes)%2==0:
        select_pairs=copy.deepcopy(selects)
    else:
        select_pairs=copy.deepcopy(selects[0])
    rank_inds=[]
    ranks=[]

    for pair in select_pairs:
        pair_ranks=dict_pair[pair.replace('.csv','')]
        ranks.append(pair_ranks)
        rank_inds.append([i for i in range(0, len(pair_ranks['select_type']))])
    string_exe1='['
    string_exe2=''
    for i in range(0,len(rank_inds)):
        if i<len(rank_inds)-1:
            string_exe1=string_exe1+'v'+str(i)+','

        else:
            string_exe1 = string_exe1 + 'v' + str(i)+']'
        string_exe2=string_exe2+'for v'+str(i)+ ' in '+'rank_inds['+str(i)+'] '
    scope = locals()
    string_exe = '[' + string_exe1 + string_exe2 + ']'

    select=eval(string_exe,scope)
    for i in select:
        type = []
        sgRNA = []
        gene_all = []
        chr = []
        strand = []
        pos = []
        on_eff = []
        off_num = []
        cut_diff = []
        for j in range(0, len(i)):
            type.append(ranks[j]['select_type'][i[j]])
            type.append(ranks[j]['select_type'][i[j]])
            if ranks[j]['select_type'][i[j]] == 'independent design':
                sgRNA.append(ranks[j]['spacer1'][i[j]])
                sgRNA.append(ranks[j]['spacer2'][i[j]])
                off_num.append(ranks[j]['off_num1'][i[j]])
                off_num.append(ranks[j]['off_num2'][i[j]])
            else:
                sgRNA.append(ranks[j]['spacer1'][i[j]])
                off_num.append(int(ranks[j]['off_num1'][i[j]]-1))
            gene_all.append(ranks[j]['gene1'][i[j]])
            gene_all.append(ranks[j]['gene2'][i[j]])
            chr.append(ranks[j]['gene1_chr'][i[j]])
            chr.append(ranks[j]['gene2_chr'][i[j]])
            strand.append(ranks[j]['gene1_strand'][i[j]])
            strand.append(ranks[j]['gene2_strand'][i[j]])
            pos.append(ranks[j]['gene1_position'][i[j]])
            pos.append(ranks[j]['gene2_position'][i[j]])
            on_eff.append(ranks[j]['cut_eff_g1'][i[j]])
            on_eff.append(ranks[j]['cut_eff_g2'][i[j]])
            cut_diff.append(ranks[j]['diff1'][i[j]])
            cut_diff.append(ranks[j]['diff2'][i[j]])
        if len(genes) % 2 > 0:
            single_gene=selects[1]
            ind=genes.index(single_gene)
            single_rank=sgss[ind]
            type.append('independent desing')
            sgRNA.append(single_rank['spacer'][0])
            gene_all.append(single_gene)
            chr.append(chr_gs[ind])
            strand.append(single_rank['strand'][0])
            pos.append(single_rank['position'][0])
            on_eff.append(float(single_rank['on_score'][0].replace('[','').replace(']','')))
            off_num.append(single_rank['num_off'][0])
            cut_diff.append(single_rank['mean_cut_diff'][0])
        types.append(type)
        sgRNAs.append(sgRNA)
        gene_alls.append(gene_all)
        chrs.append(chr)
        strands.append(strand)
        poss.append(pos)
        on_effs.append(on_eff)
        off_nums.append(off_num)
        cut_diffs.append(cut_diff)
    return types,sgRNAs,gene_alls,chrs,strands,poss,on_effs,off_nums,cut_diffs

def obtain_group_gene_info(genes,gene_all):
    two_gene_dir = rootdir+'data/two_gene_editing/'
    if len(genes) > 1:
        pairs = [genes[i] + '_' + genes[j] + '.csv' for i in range(0, len(genes) - 1) for j in range(i + 1, len(genes))]

        dict_pair = {}
        for i in range(0, len(pairs)):
            if os.path.isfile(two_gene_dir + pairs[i]):
                selects = pd.read_csv(two_gene_dir + pairs[i], header=0)
                dict_pair.update({pairs[i].replace('.csv', ''): selects})
            else:
                gene = pairs[i].replace('.csv', '').split('_')
                pair = gene[1] + '_' + gene[0] + '.csv'
                selects = pd.read_csv(two_gene_dir + pair, header=0)
                dict_pair.update({pairs[i].replace('.csv', ''): selects})
    chr_gs = []
    gene_starts = []
    gene_ends = []
    sgss = []
    for gene in genes:
        chr_g, gene_start, gene_end, sgs = gene_location(gene, gene_all)
        chr_gs.append(chr_g)
        gene_starts.append(gene_start)
        gene_ends.append(gene_ends)
        sgss.append(sgs)

    return chr_gs, gene_starts, gene_ends, sgss, dict_pair

def common_combine(commons,c_genes,c_sgRNAs,c_chr,c_strand,c_pos,c_eff,c_num,c_diff,
                   genes,select_types,select_sgRNAs,select_genes,select_cut_chrs,
                   select_cut_strands,select_cut_poss,select_cut_effs,select_off_tar_nums,
                   select_cut_diffs,type):
    for i in range(0, len(commons)):
        comm_gene = c_genes[i]
        if len(comm_gene) == len(genes):
            select_types.append([type for ot in comm_gene])
            select_sgRNAs.append(c_sgRNAs[i])
            select_genes.append(c_genes[i])
            select_cut_chrs.append(c_chr[i])
            select_cut_strands.append(c_strand[i])
            select_cut_poss.append(c_pos[i])
            select_cut_effs.append(c_eff[i])
            select_off_tar_nums.append(c_num[i])
            select_cut_diffs.append(c_diff[i])
        else:
            all_genes = copy.deepcopy(genes)
            remain_genes = list(set(all_genes).difference(comm_gene))
            select_type = [type for ot in comm_gene]
            select_sgRNA = c_sgRNAs[i]
            select_gene = c_genes[i]
            select_cut_chr = c_chr[i]
            select_cut_strand = c_strand[i]
            select_cut_pos = c_pos[i]
            select_cut_eff = c_eff[i]
            select_off_tar_num = c_num[i]
            select_cut_diff = c_diff[i]

            if len(remain_genes) == 1:
                chr_gs, gene_starts, gene_ends, sgss, dict_pair = obtain_group_gene_info(remain_genes, gene_all)
                select_type.append('independent design')
                select_sgRNA.append(sgss[0]['spacer'][0])
                select_gene.append(remain_genes[0])
                select_cut_chr.append(chr_gs[0])
                select_cut_strand.append(sgss[0]['strand'][0])
                select_cut_pos.append(sgss[0]['position'][0])
                select_cut_eff.append(sgss[0]['on_score'][0])
                select_off_tar_num.append(sgss[0]['num_off'][0])
                select_cut_diff.append(sgss[0]['mean_cut_diff'][0])
                select_types.extend(select_type)
                select_sgRNAs.extend(select_sgRNA)
                select_genes.extend(select_gene)
                select_cut_chrs.extend(select_cut_chr)
                select_cut_strands.extend(select_cut_strand)
                select_cut_poss.extend(select_cut_pos)
                select_cut_effs.extend(select_cut_eff)
                select_off_tar_nums.extend(select_off_tar_num)
                select_cut_diffs.extend(select_cut_diff)
            else:
                combs_selects = select_combination_all(remain_genes)
                chr_gs, gene_starts, gene_ends, sgss, dict_pair = obtain_group_gene_info(remain_genes, gene_all)
                for j in range(0, len(combs_selects)):
                    types, sgRNAs, gene_alls, chrs, strands, poss, on_effs, off_nums, cut_diffs = comine_pairs(sgss,
                                                                                                               dict_pair,
                                                                                                               combs_selects[
                                                                                                                   i],
                                                                                                               remain_genes,
                                                                                                               chr_gs)
                    for k in range(0, len(sgRNAs)):
                        types[k].extend(select_type)
                        sgRNAs[k].extend(select_sgRNA)
                        gene_alls[k].extend(select_gene)
                        chrs[k].extend(select_cut_chr)
                        strands[k].extend(select_cut_strand)
                        poss[k].extend(select_cut_pos)
                        on_effs[k].extend(select_cut_eff)
                        off_nums[k].extend(select_off_tar_num)
                        cut_diffs[k].extend(select_cut_diff)
                    select_types.extend(types)
                    select_sgRNAs.extend(sgRNAs)
                    select_genes.extend(gene_alls)
                    select_cut_chrs.extend(chrs)
                    select_cut_strands.extend(strands)
                    select_cut_poss.extend(poss)
                    select_cut_effs.extend(on_effs)
                    select_off_tar_nums.extend(off_nums)
                    select_cut_diffs.extend(cut_diffs)
    return select_types,select_sgRNAs,select_genes,select_cut_chrs,select_cut_strands,select_cut_poss,select_cut_effs,select_off_tar_nums,select_cut_diffs

def multi_gene(input_genes, gene_all):
    genes=input_genes.strip().split(',')
    comm_genes, comm_sgRNAs, comm_cut_eff, comm_off_num, comm_cut_diff, comm_cut_strand, comm_cut_chr, comm_cut_pos=common_one2all(genes,gene_all)
    comb_genes, comb_sgRNAs, comb_cut_eff, comb_off_num, comb_cut_diff, comb_cut_strand, comb_cut_chr, comb_cut_pos=off_one2all(genes)
    select_types=[]
    select_sgRNAs = []
    select_genes = []
    select_cut_chrs = []
    select_cut_strands = []
    select_cut_poss = []
    select_cut_effs = []
    select_off_tar_nums = []
    select_cut_diffs = []
    if len(comm_sgRNAs)+len(comb_sgRNAs)==0:
        combs_selects = select_combination_all(genes)
        chr_gs, gene_starts, gene_ends, sgss, dict_pair=obtain_group_gene_info(genes,gene_all)
        for i in range(0,len(combs_selects)):
            types, sgRNAs, gene_alls, chrs, strands, poss, on_effs, off_nums, cut_diffs=comine_pairs(sgss, dict_pair, combs_selects[i], genes, chr_gs)
            select_types.extend(types)
            select_sgRNAs.extend(sgRNAs)
            select_genes.extend(gene_alls)
            select_cut_chrs.extend(chrs)
            select_cut_strands.extend(strands)
            select_cut_poss.extend(poss)
            select_cut_effs.extend(on_effs)
            select_off_tar_nums.extend(off_nums)
            select_cut_diffs.extend(cut_diffs)
    elif len(comm_sgRNAs)>0 and len(comb_sgRNAs)==0:
        for i in range(0,len(comm_sgRNAs)):
            comm_gene=comm_genes[i]
            if len(comm_gene)==len(genes):
                select_types.append(['common' for ot in comm_gene])
                select_sgRNAs.append(comm_sgRNAs[i])
                select_genes.append(comm_genes[i])
                select_cut_chrs.append(comm_cut_chr[i])
                select_cut_strands.append(comm_cut_strand[i])
                select_cut_poss.append(comm_cut_pos[i])
                select_cut_effs.append(comm_cut_eff[i])
                select_off_tar_nums.append(comm_off_num[i])
                select_cut_diffs.append(comm_cut_diff[i])
            else:
                all_genes=copy.deepcopy(genes)
                remain_genes=list(set(all_genes).difference(comm_gene))
                select_type=['common' for ot in comm_gene]
                select_sgRNA=comm_sgRNAs[i]
                select_gene = comm_genes[i]
                select_cut_chr = comm_cut_chr[i]
                select_cut_strand = comm_cut_strand[i]
                select_cut_pos = comm_cut_pos[i]
                select_cut_eff = comm_cut_eff[i]
                select_off_tar_num = comm_off_num[i]
                select_cut_diff = comm_cut_diff[i]

                if len(remain_genes)==1:
                    chr_gs, gene_starts, gene_ends, sgss, dict_pair = obtain_group_gene_info(remain_genes, gene_all)
                    select_type.append('independent design')
                    select_sgRNA.append(sgss[0]['spacer'][0])
                    select_gene.append(remain_genes[0])
                    select_cut_chr.append(chr_gs[0])
                    select_cut_strand.append(sgss[0]['strand'][0])
                    select_cut_pos.append(sgss[0]['position'][0])
                    select_cut_eff.append(sgss[0]['on_score'][0])
                    select_off_tar_num.append(sgss[0]['num_off'][0])
                    select_cut_diff.append(sgss[0]['mean_cut_diff'][0])
                    select_types.extend(select_type)
                    select_sgRNAs.extend(select_sgRNA)
                    select_genes.extend(select_gene)
                    select_cut_chrs.extend(select_cut_chr)
                    select_cut_strands.extend(select_cut_strand)
                    select_cut_poss.extend(select_cut_pos)
                    select_cut_effs.extend(select_cut_eff)
                    select_off_tar_nums.extend(select_off_tar_num)
                    select_cut_diffs.extend(select_cut_diff)
                else:
                    combs_selects = select_combination_all(remain_genes)
                    chr_gs, gene_starts, gene_ends, sgss, dict_pair = obtain_group_gene_info(remain_genes, gene_all)
                    for j in range(0, len(combs_selects)):
                        types, sgRNAs, gene_alls, chrs, strands, poss, on_effs, off_nums, cut_diffs = comine_pairs(sgss, dict_pair, combs_selects[i], remain_genes, chr_gs)
                        for k in range(0,len(sgRNAs)):
                            types[k].extend(select_type)
                            sgRNAs[k].extend(select_sgRNA)
                            gene_alls[k].extend(select_gene)
                            chrs[k].extend(select_cut_chr)
                            strands[k].extend(select_cut_strand)
                            poss[k].extend(select_cut_pos)
                            on_effs[k].extend(select_cut_eff)
                            off_nums[k].extend(select_off_tar_num)
                            cut_diffs[k].extend(select_cut_diff)
                        select_types.extend(types)
                        select_sgRNAs.extend(sgRNAs)
                        select_genes.extend(gene_alls)
                        select_cut_chrs.extend(chrs)
                        select_cut_strands.extend(strands)
                        select_cut_poss.extend(poss)
                        select_cut_effs.extend(on_effs)
                        select_off_tar_nums.extend(off_nums)
                        select_cut_diffs.extend(cut_diffs)
    elif len(comm_sgRNAs)==0 and len(comb_sgRNAs)>0:
        c_sgRNAs = comb_sgRNAs
        c_genes = comb_genes
        c_chr = comb_cut_chr
        c_strand = comb_cut_strand
        c_pos = comb_cut_pos
        c_eff = comb_cut_eff
        c_num = comb_off_num
        c_diff = comb_cut_diff
        type = 'off target'
        select_types, select_sgRNAs, select_genes, select_cut_chrs, \
        select_cut_strands, select_cut_poss, select_cut_effs, \
        select_off_tar_nums, select_cut_diffs = common_combine(c_sgRNAs, c_genes, c_sgRNAs, c_chr, c_strand, c_pos, c_eff, c_num,
                                                               c_diff, genes,
                                                               select_types, select_sgRNAs, select_genes,
                                                               select_cut_chrs, select_cut_strands, select_cut_poss,
                                                               select_cut_effs,
                                                               select_off_tar_nums, select_cut_diffs, type)
    elif len(comm_sgRNAs)>0 and len(comb_sgRNAs)>0:
        for i in range(0, len(comm_sgRNAs)):
            for j in range(0, len(comb_sgRNAs)):
                intersect_com=np.intersect1d(comm_sgRNAs[i],comb_sgRNAs[j])
                if len(intersect_com)>0 :
                    if comm_sgRNAs[i]==comb_sgRNAs[j] or len(comm_sgRNAs[i])>=len(comb_sgRNAs[j]):
                        c_sgRNAs=[comm_sgRNAs[i]]
                        c_genes=[comm_genes[i]]
                        c_chr=[comm_cut_chr[i]]
                        c_strand=[comm_cut_strand[i]]
                        c_pos=[comm_cut_pos[i]]
                        c_eff=[comm_cut_eff[i]]
                        c_num=[comm_off_num[i]]
                        c_diff=[comm_cut_diff[i]]
                        type='common'

                    else:
                        c_sgRNAs = [comb_sgRNAs[j]]
                        c_genes = [comb_genes[j]]
                        c_chr = [comb_cut_chr[j]]
                        c_strand = [comb_cut_strand[j]]
                        c_pos = [comb_cut_pos[j]]
                        c_eff = [comb_cut_eff[j]]
                        c_num = [comb_off_num[j]]
                        c_diff = [comb_cut_diff[j]]
                        type='off target'

                    select_types, select_sgRNAs, select_genes, select_cut_chrs, \
                        select_cut_strands, select_cut_poss, select_cut_effs, \
                        select_off_tar_nums, select_cut_diffs=common_combine(c_sgRNAs, c_genes,c_sgRNAs,c_chr,c_strand,c_pos,c_eff,c_num,c_diff, genes,
                                                                             select_types,select_sgRNAs,select_genes,
                                                                            select_cut_chrs, select_cut_strands,select_cut_poss,select_cut_effs,
                                                                            select_off_tar_nums, select_cut_diffs,type)
                else:
                    select_type=['common' for t in comm_genes[i]]
                    select_sgRNA=copy.deepcopy(comm_sgRNAs[i])
                    select_gene=copy.deepcopy(comm_genes[i])
                    select_cut_chr=copy.deepcopy(comm_cut_chr[i])
                    select_cut_strand=copy.deepcopy(comm_cut_strand[i])
                    select_cut_pos=copy.deepcopy(comm_cut_pos[i])
                    select_cut_eff=copy.deepcopy(comm_cut_eff[i])
                    select_off_tar_num=copy.deepcopy(comm_off_num[i])
                    select_cut_diff=copy.deepcopy(comm_cut_diff[i])
                    select_type.extend(['off target' for t in comb_genes[j]])
                    select_sgRNA.extend(comb_sgRNAs[j])
                    select_gene.extend(comb_genes[j])
                    select_cut_chr.extend(comb_cut_chr[j])
                    select_cut_strand.extend(comb_cut_strand[j])
                    select_cut_pos.extend(comb_cut_pos[j])
                    select_cut_eff.extend(comb_cut_eff[j])
                    select_off_tar_num.extend(comm_off_num[j])
                    select_cut_diff.extend(comm_cut_diff[j])

                    all_genes = copy.deepcopy(genes)
                    remain_genes = list(set(all_genes).difference(select_gene))

                    if len(remain_genes) == 1:
                        chr_gs, gene_starts, gene_ends, sgss, dict_pair = obtain_group_gene_info(remain_genes, gene_all)
                        select_type.append('independent design')
                        select_sgRNA.append(sgss[0]['spacer'][0])
                        select_gene.append(remain_genes[0])
                        select_cut_chr.append(chr_gs[0])
                        select_cut_strand.append(sgss[0]['strand'][0])
                        select_cut_pos.append(sgss[0]['position'][0])
                        select_cut_eff.append(sgss[0]['on_score'][0])
                        select_off_tar_num.append(sgss[0]['num_off'][0])
                        select_cut_diff.append(sgss[0]['mean_cut_diff'][0])
                        select_types.extend(select_type)
                        select_sgRNAs.extend(select_sgRNA)
                        select_genes.extend(select_gene)
                        select_cut_chrs.extend(select_cut_chr)
                        select_cut_strands.extend(select_cut_strand)
                        select_cut_poss.extend(select_cut_pos)
                        select_cut_effs.extend(select_cut_eff)
                        select_off_tar_nums.extend(select_off_tar_num)
                        select_cut_diffs.extend(select_cut_diff)
                    else:
                        combs_selects = select_combination_all(remain_genes)
                        chr_gs, gene_starts, gene_ends, sgss, dict_pair = obtain_group_gene_info(remain_genes, gene_all)
                        for j in range(0, len(combs_selects)):
                            types, sgRNAs, gene_alls, chrs, strands, poss, on_effs, off_nums, cut_diffs = comine_pairs(
                                sgss, dict_pair, combs_selects[i], remain_genes, chr_gs)
                            for k in range(0, len(sgRNAs)):
                                types[k].extend(select_type)
                                sgRNAs[k].extend(select_sgRNA)
                                gene_alls[k].extend(select_gene)
                                chrs[k].extend(select_cut_chr)
                                strands[k].extend(select_cut_strand)
                                poss[k].extend(select_cut_pos)
                                on_effs[k].extend(select_cut_eff)
                                off_nums[k].extend(select_off_tar_num)
                                cut_diffs[k].extend(select_cut_diff)
                            select_types.extend(types)
                            select_sgRNAs.extend(sgRNAs)
                            select_genes.extend(gene_alls)
                            select_cut_chrs.extend(chrs)
                            select_cut_strands.extend(strands)
                            select_cut_poss.extend(poss)
                            select_cut_effs.extend(on_effs)
                            select_off_tar_nums.extend(off_nums)
                            select_cut_diffs.extend(cut_diffs)

    chr_gs, gene_starts, gene_ends, sgss, dict_pair = obtain_group_gene_info(genes, gene_all)
    select_types.append(['independent design' for s in genes])
    select_sgRNAs.append([sgss[s]['spacer'][0] for s in range(0, len(genes))])
    select_genes.append(genes)
    select_cut_chrs.append(chr_gs)
    select_cut_strands.append([sgss[s]['strand'][0] for s in range(0, len(genes))])
    select_cut_poss.append([sgss[s]['position'][0] for s in range(0, len(genes))])
    select_cut_effs.append([float(sgss[s]['on_score'][0].replace('[','').replace(']','')) for s in range(0, len(genes))])
    select_off_tar_nums.append([sgss[s]['num_off'][0] for s in range(0, len(genes))])
    select_cut_diffs.append([sgss[s]['mean_cut_diff'][0] for s in range(0, len(genes))])
    return select_types,select_sgRNAs,select_genes,select_cut_chrs,select_cut_strands,select_cut_poss,select_cut_effs,select_off_tar_nums,select_cut_diffs

def multi_gene_rank_selects(input_genes, select_types,select_sgRNAs,select_genes,select_cut_chrs,select_cut_strands,select_cut_poss,select_cut_effs,select_off_tar_nums,select_cut_diffs):
    effs=[]
    off_nums=[]
    diffs=[]
    for i in range(0, len(select_sgRNAs)):
        effs.append(statistics.harmonic_mean(select_cut_effs[i]))
        off_nums.append(sum(select_off_tar_nums[i]))
        diff=[]
        s=select_cut_diffs[i]
        for j in s:
            if j<0:
                diff.append(0.0000001)
            else:
                diff.append(j)
        diffs.append(statistics.harmonic_mean(diff))

    rank_avg_cut = ranks(effs, True)
    rank_off = ranks(off_nums, False)
    rank_diff = ranks(diffs, True)
    final_rank = []
    for s in range(0, len(effs)):
        rank = (rank_avg_cut[s][2] + rank_off[s][2] + rank_diff[s][2]) / 3
        final_rank.append(rank)

    final_index = np.zeros((len(final_rank), 2))
    for i in range(0, len(final_index)):
        final_index[i, 0] = final_rank[i]
        final_index[i, 1] = i
    # Out = np.hstack((final_rank, final_index))
    out = sorted(final_index, key=lambda row: row[0], reverse=False)
    file_path = rootdir+'data/multi_genes/' + input_genes+ '.csv'
    with open(file_path, 'w') as csv_file:
        fieldnames = ['select_type','select_spacers', 'cut_genes', 'gene_chrs', 'gene_strands',
                      'gene_positions', 'cut_effs', 'average_cut_on',
                      'rank_on', 'off_nums', 'total_off-tar_num', 'rank_off', 'diffs', 'mean_diff', 'rank_diff']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(0, len(effs)):
            ind = int(out[i][1])
            writer.writerow(
                {'select_type':select_types[ind],'select_spacers': select_sgRNAs[ind],
                 'cut_genes':  select_genes[ind], 'gene_chrs': select_cut_chrs[ind],
                 'gene_strands': select_cut_strands[ind], 'gene_positions': select_cut_poss[ind],
                 'cut_effs': select_cut_effs[ind],
                 'average_cut_on': effs[ind], 'rank_on': rank_avg_cut[ind][2],'off_nums':select_off_tar_nums[ind],
                 'total_off-tar_num': off_nums[ind], "rank_off": rank_off[ind][2],
                 'diffs': select_cut_diffs[ind], 'mean_diff': diffs[ind],
                 'rank_diff': rank_diff[ind][2]})
    return file_path

def common_one2all(genes,gene_all):
    spacers_all=[]
    chr_gs=[]
    sgss=[]
    comm_genes = []
    comm_sgRNAs = []
    comm_cut_eff = []
    comm_off_num = []
    comm_cut_diff = []
    comm_cut_strand = []
    comm_cut_chr = []
    comm_cut_pos = []
    for gene in genes:
        sgRNAs_ex, spacers, strand, poss, pre_on_scores=readOnTars(gene)
        chr_g, gene_start, gene_end, sgs = gene_location(gene, gene_all)
        spacers_all.append(spacers)
        sgss.append(sgs)
        chr_gs.append(chr_g)

    n=len(genes)

    index=[j for j in range(0,n)]
    combs=[]
    commons=[]
    for i in range(3,n+1):
        comb=list(iter.combinations(index, i))

        for j in comb:
            ind=[k for k in j]
            comm_gene = [genes[k] for k in j]
            common=reduce(np.intersect1d, ([spacers_all[k] for k in j]))
            cut_gene = []
            cut_sgRNA = []
            cut_chr = []
            cut_strand = []
            cut_pos = []
            cut_eff = []
            cut_off_num = []
            cut_diff = []
            if len(common)>0:
                for t in common:
                    flag=0
                    for m in j:
                        ind_ori = spacers_all[m].index(t)

                        cut_gene.append(genes[m])
                        if flag == 0:
                            cut_sgRNA.append(t)
                            cut_off_num.append(sgss[m]['num_off'][ind_ori])
                            flag = flag + 1
                        cut_chr.append(chr_gs[m])
                        cut_strand.append(sgss[m]['strand'][ind_ori])
                        cut_pos.append(sgss[m]['position'][ind_ori])
                        cut_eff.append(float(sgss[m]['on_score'][ind_ori].replace('[','').replace(']','')))
                        cut_diff.append(sgss[m]['mean_cut_diff'][ind_ori])
                comm_genes.append(cut_gene)
                comm_sgRNAs.append(cut_sgRNA)
                comm_cut_eff.append(cut_eff)
                comm_off_num.append(cut_off_num)
                comm_cut_diff.append(cut_diff)
                comm_cut_strand.append(cut_strand)
                comm_cut_chr.append(cut_chr)
                comm_cut_pos.append(cut_pos)
    return comm_genes, comm_sgRNAs, comm_cut_eff, comm_off_num, comm_cut_diff, comm_cut_strand, comm_cut_chr, comm_cut_pos


def off_one2all(genes):
    comb_genes=[]
    comb_sgRNAs=[]
    comb_cut_eff=[]
    comb_off_num=[]
    comb_cut_diff=[]
    comb_cut_strand=[]
    comb_cut_chr=[]
    comb_cut_pos=[]

    off_anno_dir=rootdir+'data/single_gene_MPCP1/'
    for gene in genes:
        off_anno=pd.read_csv(off_anno_dir+gene+'_0.5_MPCP.csv', header=0)
        off_gene=[]
        off_spacers=[]
        tar_gene=[]
        tar_spacer=[]
        #ontar_genes.append(gene)
        tar_index=[]
        rank_ind=[]
        for i in range(0,len(off_anno['gene2'])):
            #print(str(off_anno['gene2']))
            if str(off_anno['gene2'][i])=='nan':
                continue
            else:
                ges=str(off_anno['gene2'][i]).strip().split('|')
                off_gene.extend(ges)
                off_spacer=[off_anno['select_spacer'][i] for j in ges]
                off_spacers.extend(off_spacer)
                rank_ind.extend([i for j in ges])
        gene_remain=copy.deepcopy(genes)
        gene_remain.remove(gene)
        for remain in gene_remain:
            index=np.array([i for i, X in enumerate(off_gene) if X == remain])
            if len(index)>0:
               tar_gene.append(remain)
               tar_spacer.append([off_spacers[k] for k in index])
               tar_index.append([rank_ind[k] for k in index])
        if len(tar_spacer)>0:
            inds=[k for k in range(0, len(tar_gene))]
            for i in range(2, len(tar_gene) + 1):
                comb = list(iter.combinations(inds, i))

                for j in comb:
                    comb_gene = [tar_gene[k] for k in j]
                    common = reduce(np.intersect1d, ([tar_spacer[k] for k in j]))
                    cut_gene=[]
                    cut_sgRNA=[]
                    cut_chr=[]
                    cut_strand=[]
                    cut_pos=[]
                    cut_eff=[]
                    cut_off_num=[]
                    cut_diff=[]
                    if len(common) > 0:
                        for t in common:
                            cut_gene.append(gene)
                            cut_gene.extend(comb_gene)
                            cut_sgRNA.append(t)
                            cut_chr.append(off_anno['gene1_chr'][0])
                            cut_strand.append(off_anno['gene1_strand'][0])
                            cut_pos.append(off_anno['gene1_position'][0])
                            flag=0
                            for m in j:
                                ind=tar_spacer[m].index(t)
                                ind_ori=tar_index[m][ind]

                                cut_chr.append(off_anno['gene2_chr'][ind_ori])
                                cut_strand.append(off_anno['gene2_strand'][ind_ori])
                                cut_pos.append(off_anno['gene2_position'][ind_ori])
                                if flag==0:
                                    cut_eff.append(off_anno['cut_eff_g1'][ind_ori])
                                    cut_diff.append(off_anno['on_diff'][ind_ori])
                                    cut_off_num.append(int(off_anno['total_off-tar_num'][ind_ori]) - len(comb_gene))
                                    flag=flag+1
                                cut_eff.append(off_anno['cut_eff_g2'][ind_ori])

                                cut_diff.append(off_anno['off_diff'][ind_ori])
                        comb_genes.append(cut_gene)
                        comb_sgRNAs.append(cut_sgRNA)
                        comb_cut_eff.append(cut_eff)
                        comb_off_num.append(cut_off_num)
                        comb_cut_diff.append(cut_diff)
                        comb_cut_strand.append(cut_strand)
                        comb_cut_chr.append(cut_chr)
                        comb_cut_pos.append(cut_pos)
    return comb_genes,comb_sgRNAs,comb_cut_eff,comb_off_num,comb_cut_diff,comb_cut_strand,comb_cut_chr,comb_cut_pos


def select_combination_all(genes):
    #genes = input_genes.strip().split(',')
    n_g=len(genes)
    if (n_g % 2)==0:
        pairs=[genes[i]+'_'+genes[j]+'.csv' for i in range(0,len(genes)-1) for j in range(i+1,len(genes))]
        uni_genes=copy.deepcopy(genes)
        combs_selects=select_comb_even(pairs,uni_genes)
    else:
        combs_selects=[]
        for gene in genes:
            genes0=copy.deepcopy(genes)
            genes0.remove(gene)
            #genes0=genes0.remove(gene)
            pairs=[genes0[i]+'_'+genes0[j]+'.csv' for i in range(0,len(genes0)-1) for j in range(i+1,len(genes0))]
            combs_select = select_comb_even(pairs, genes0)
            for coms in combs_select:
                com=[]
                com.append(coms)
                com.append(gene)
                combs_selects.append(com)
    return combs_selects

def pair_gene(pair):
    gene=pair.replace('.csv','').split('_')
    return gene

def select_comb_even(pairs,uni_genes):
    n_g=len(uni_genes)
    #n_p=len(pairs)
    combs_select=[]
    select_num=int(n_g/2)
    #select_genes=copy.deepcopy(uni_genes)
    combs=list(iter.combinations(pairs, select_num))
    for i in combs:
        gene_comb=[]
        for j in range(0,len(i)):
            gene=pair_gene(i[j])
            gene_comb.extend(gene)
        gene_uni=np.unique(gene_comb)
        if len(gene_uni)==len(gene_comb):
            combs_select.append(i)
    return combs_select

def multi_gene_selection(input_genes):
    anno_file = rootdir+'data/Homo_sapiens.GRCh38.92.gtf'
    gene_chr, gene_start, gene_end, gene_strand, gene_id, gene_name = obtain_exons(anno_file)
    gene_all = pd.DataFrame(zip(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand),
                            columns=['id', 'name', 'chr', 'start', 'end', 'strand'])
    ids, symbols, chrs, starts, ends, strands = obtain_gene_info()
    select_types, select_sgRNAs, select_genes, select_cut_chrs, select_cut_strands, select_cut_poss, select_cut_effs, select_off_tar_nums, select_cut_diffs = multi_gene(
        input_genes, gene_all)
    file_path = multi_gene_rank_selects(input_genes, select_types, select_sgRNAs, select_genes, select_cut_chrs,
                                        select_cut_strands,
                                        select_cut_poss, select_cut_effs, select_off_tar_nums, select_cut_diffs)
    return file_path

if __name__ == '__main__':
    #name = 'cancer_genes'
    rootdir = '../'
    anno_file = rootdir+'data/Homo_sapiens.GRCh38.92.gtf'
    gene_chr, gene_start, gene_end, gene_strand, gene_id, gene_name = obtain_exons(anno_file)
    gene_all = pd.DataFrame(zip(gene_id, gene_name, gene_chr, gene_start, gene_end, gene_strand),
                            columns=['id', 'name', 'chr', 'start', 'end', 'strand'])
    ids, symbols, chrs, starts, ends, strands = obtain_gene_info()

    #input_genes='ENSG00000148584,ENSG00000097007,ENSG00000164398,ENSG00000009709,ENSG00000005073,ENSG00000128713'
    #input_genes = 'ENSG00000148584,ENSG00000097007,ENSG00000164398'
    input_genes = 'ENSG00000148584,ENSG00000097007,ENSG00000164398,ENSG00000009709,ENSG00000005073'
    #input_genes = 'ENSG00000148584,ENSG00000097007,ENSG00000164398,ENSG00000009709'
    select_types, select_sgRNAs, select_genes, select_cut_chrs, select_cut_strands, select_cut_poss, select_cut_effs, select_off_tar_nums, select_cut_diffs=multi_gene(input_genes, gene_all)
    file_path=multi_gene_rank_selects(input_genes, select_types, select_sgRNAs, select_genes, select_cut_chrs, select_cut_strands,
                            select_cut_poss, select_cut_effs, select_off_tar_nums, select_cut_diffs)




