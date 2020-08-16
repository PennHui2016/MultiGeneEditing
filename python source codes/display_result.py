import sys
#sys.path.append("D:/project/quantification/code/GAN_learn/venv/Lib/site-packages")
import pandas as pd
import numpy as np
import multi_gene_editing_offline as mg
import csv
import time
import os


#rootdir=os.path.abspath('..').replace('\\','/')+'/'
#single_gene_dir=rootdir+'data/single_editing_modify/'
#single_gene_mpcp_dir=rootdir+'data/single_gene_MPCP1/'
#two_gene_dir=rootdir+'data/two_gene_editing/'
#multi_gene_dir=rootdir+'data/multi_genes/'

def single_gene_mpcp(gene,Time):
    file=single_gene_mpcp_dir+gene+'_0.5_MPCP.csv'
    records=pd.read_csv(file)
    mpcp_genes=records['gene2']

    ind=0
    for i in range(0,len(mpcp_genes)):
        s=str(mpcp_genes[i])
        if s=='nan':
            continue
        elif s==gene:
            continue
        else:
            ind=i
            break

    select_MPCP=records[ind:ind+1]
    #Time=time.time()
    file_path = rootdir+'results/single_gene_MPCP_'+str(Time)+'.csv'
    with open(file_path, 'w',newline='') as csv_file:
        fieldnames = ['Spacer', 'Spacer ex', 'strand', 'position', 'partner',
                      'p_chr', 'p_pos', 'p_stra',
                      'partner site', 'cut_on', 'cut_partner', 'total off num', 'on_PCS', 'partner_PCS', 'mean_PCS']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(
            {'Spacer': select_MPCP['select_spacer'][ind].replace(',',';'), 'Spacer ex': select_MPCP['select_spacer_ex'][ind].replace(',',';'), 'strand': select_MPCP['gene1_strand'][ind].replace(',',';'),
             'position': str(select_MPCP['gene1_position'][ind]).replace(',',';'), 'partner': select_MPCP['gene2'][ind].replace(',',';'),
             'p_chr': str(select_MPCP['gene2_chr'][ind]).replace(',',';'), 'p_pos': str(select_MPCP['gene2_position'][ind]).replace(',',';'),
             'p_stra': select_MPCP['gene2_strand'][ind].replace(',',';'), 'partner site': select_MPCP['off_seq'][ind].replace(',',';'),
             'cut_on': str(select_MPCP['cut_eff_g1'][ind]).replace(',',';'), 'cut_partner': str(select_MPCP['cut_eff_g2'][ind]).replace(',',';'), 'total off num': str(select_MPCP['total_off-tar_num'][ind]).replace(',',';'),
             'on_PCS': str(select_MPCP['on_diff'][ind]).replace(',',';'), 'partner_PCS': str(select_MPCP['off_diff'][ind]).replace(',',';'), 'mean_PCS': str(select_MPCP['mean_diff'][ind]).replace(',',';')})
    return file_path,select_MPCP

def single_gene_select_spacer(gene,Time):

    file=single_gene_dir+gene+'_0.5.csv'
    records = pd.read_csv(file)

    n=len(records['spacer'])
    if n<5:
        N=n
    else:
        N=5

    selected_spacers=records[0:N]
    #Time=time.time()
    file_path = rootdir+'results/single_gene_spacer_'+str(Time)+'.csv'
    with open(file_path, 'w',newline='') as csv_file:
        fieldnames = ['index', 'Spacer sequence', 'Spacer extend sequence', 'position', 'strand',
                      'cut efficiency', 'off target site number', 'mean PCS']

        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(0, N):
            writer.writerow(
            {'index': i, 'Spacer sequence': selected_spacers['spacer'][i].replace(',',';'),
             'Spacer extend sequence': selected_spacers['spacer_ex'][i].replace(',',';'),
             'position': str(selected_spacers['position'][i]).replace(',',';'), 'strand': str(selected_spacers['strand'][i]).replace(',',';'),
             'cut efficiency': str(selected_spacers['on_score'][i]).replace(',',';'), 'off target site number': str(selected_spacers['num_off'][i]).replace(',',';'),
             'mean PCS': str(selected_spacers['mean_cut_diff'][i]).replace(',',';')})
    return file_path,selected_spacers

def two_gene_select_spacer(input_gene,Time):
    input_gene=input_gene.replace(',','_')
    file = two_gene_dir + input_gene + '.csv'
    records = pd.read_csv(file)
    select_types=list(records['select_type'])

    select_inds = []
    select_Types=[]
    spacer_num=[]
    spacers=[]
    genes=[]
    chrs=[]
    strands=[]
    positions=[]
    cut_ons=[]
    diffs=[]

    if 'independent design' in select_types:
        ind_indep = select_types.index('independent design')
        select_inds.append(ind_indep)
        select_Types.append('independent design;independent design')
        spacer_num.append(2)
        sp=records['spacer1'][ind_indep]+';'+records['spacer2'][ind_indep]
        spacers.append(sp)
        genes.append(records['gene1'][ind_indep]+','+records['gene2'][ind_indep])
        chrs.append(str(records['gene1_chr'][ind_indep])+';'+str(records['gene2_chr'][ind_indep]))
        positions.append(str(records['gene1_position'][ind_indep])+';' + str(records['gene2_position'][ind_indep]))
        strands.append(records['gene2_strand'][ind_indep]+';' + records['gene2_strand'][ind_indep])
        cut_ons.append(str(records['cut_eff_g1'][ind_indep])+';' + str(records['cut_eff_g2'][ind_indep]))
        diffs.append(str(records['diff1'][ind_indep])+';' + str(records['diff2'][ind_indep]))

    if 'off target' in select_types:
        ind_off = select_types.index('off target')
        select_inds.append(ind_off)
        select_Types.append('off target;off target')
        spacer_num.append(1)
        sp = records['spacer1'][ind_off]
        spacers.append(sp)
        genes.append(records['gene1'][ind_off] + ',' + records['gene2'][ind_off])
        chrs.append(str(records['gene1_chr'][ind_off]) + ';' + str(records['gene2_chr'][ind_off]))
        positions.append(str(records['gene1_position'][ind_off]) + ';' + str(records['gene2_position'][ind_off]))
        strands.append(records['gene2_strand'][ind_off] + ';' + records['gene2_strand'][ind_off])
        cut_ons.append(str(records['cut_eff_g1'][ind_off]) + ';' + str(records['cut_eff_g2'][ind_off]))
        diffs.append(str(records['diff1'][ind_off]) + ';' + str(records['diff2'][ind_off]))

    if 'common spacer' in select_types:
        ind_comm = select_types.index('common spacer')
        select_inds.append(ind_comm)
        select_Types.append('common;common')
        spacer_num.append(1)
        sp = records['spacer1'][ind_comm]
        spacers.append(sp)
        genes.append(records['gene1'][ind_comm] + ',' + records['gene2'][ind_comm])
        chrs.append(str(records['gene1_chr'][ind_comm]) + ';' + str(records['gene2_chr'][ind_comm]))
        positions.append(str(records['gene1_position'][ind_comm]) + ';' + str(records['gene2_position'][ind_comm]))
        strands.append(records['gene2_strand'][ind_comm] + ';' + records['gene2_strand'][ind_comm])
        cut_ons.append(str(records['cut_eff_g1'][ind_comm]) + ';' + str(records['cut_eff_g2'][ind_comm]))
        diffs.append(str(records['diff1'][ind_comm]) + ';' + str(records['diff2'][ind_comm]))

    if len(select_inds)>1:
        select_spacers=records[select_inds[0]:select_inds[0]+1]
        for i in range(1,len(select_inds)):
            select_spacers=pd.concat([select_spacers, records[select_inds[i]:select_inds[i]+1]],ignore_index=True)
    elif len(select_inds)==1:
        select_spacers=records[select_inds[0]:select_inds[0]+1]

    #Time=time.time()
    file_path = rootdir+'results/multi_gene_spacer_'+str(Time)+'.csv'
    with open(file_path, 'w',newline='') as csv_file:
        fieldnames = ['spacer num', 'select_types', 'Spacers', 'genes', 'chrs',
                      'strands', 'positions', 'cut_ons', 'total off num', 'PCSs', 'average_PCS']

        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        N=len(select_spacers['select_type'])
        for i in range(0, N):
            writer.writerow(
                {'spacer num': str(spacer_num[i]).replace(',',';'), 'select_types': select_Types[i].replace(',',';'),
                 'Spacers': spacers[i].replace(',',';'),
                 'genes': genes[i].replace(',',';'), 'chrs': chrs[i].replace(',',';'),
                 'strands': strands[i].replace(',',';'),
                 'positions': positions[i].replace(',',';'),
                 'cut_ons': str(cut_ons[i]).replace(',',';'),'total off num':str(select_spacers['total_off-tar_num'][i]).replace(',',';'),'PCSs':str(diffs[i]).replace(',',';'),'average_PCS':str(select_spacers['mean_diff'][i]).replace(',',';')})
    return file_path,select_spacers

def multi_gene_select_spacer(input_genes,Time):
    file_path=mg.multi_gene_selection(input_genes)
    records=pd.read_csv(file_path)
    os.remove(file_path)
    spacer_nums=[]
    spacers=list(records['select_spacers'])
    for genes in spacers:
        ges=genes.split(',')
        spacer_nums.append(len(ges))
    uni_num=np.unique(spacer_nums)

    sp_nums=[]
    if len(uni_num)==1:
        ind=spacer_nums.index(uni_num[0])
        select_spacers=records[ind:ind+1]
        sp_nums.append(spacer_nums[ind])
    elif len(uni_num)>1:
        ind0=spacer_nums.index(uni_num[0])
        select_spacers=records[ind0:ind0+1]
        sp_nums.append(spacer_nums[ind0])
        for i in range(1,len(uni_num)):
            ind=spacer_nums.index(uni_num[i])
            select_spacers = pd.concat([select_spacers, records[ind:ind+1]], ignore_index=True)
            sp_nums.append(spacer_nums[ind])

    #Time=time.time()
    file_path = rootdir+'results/multi_gene_spacer_'+str(Time)+'.csv'
    with open(file_path, 'w',newline='') as csv_file:
        fieldnames = ['spacer num', 'select_types', 'Spacers', 'genes', 'chrs',
                      'strands', 'positions', 'cut_ons', 'total off num', 'PCSs', 'average_PCS']

        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        N = len(select_spacers['select_type'])
        for i in range(0, N):
            writer.writerow(
                {'spacer num': str(sp_nums[i]).replace(',',';'), 'select_types': select_spacers['select_type'][i].replace(',',';'),
                 'Spacers': select_spacers['select_spacers'][i].replace(',',';'),
                 'genes': select_spacers['cut_genes'][i].replace(',',';'), 'chrs': select_spacers['gene_chrs'][i].replace(',',';'),
                 'strands': select_spacers['gene_strands'][i].replace(',',';'),
                 'positions': select_spacers['gene_positions'][i].replace(',',';'),
                 'cut_ons': str(select_spacers['cut_effs'][i]).replace(',',';'), 'total off num': str(select_spacers['total_off-tar_num'][i]).replace(',',';'), 'PCSs': select_spacers['diffs'][i].replace(',',';'),
                 'average_PCS': str(select_spacers['mean_diff'][i]).replace(',',';')})
    return file_path

if __name__ == '__main__':
    input_genes=sys.argv[1]
    Time=sys.argv[2]
    rootdir=sys.argv[3].replace('\\','/')+'/'
    single_gene_dir=rootdir+'data/single_editing_modify/'
    single_gene_mpcp_dir=rootdir+'data/single_gene_MPCP1/'
    two_gene_dir=rootdir+'data/two_gene_editing/'
    multi_gene_dir=rootdir+'data/multi_genes/'
    genes=input_genes.split(',')
    gene_num=len(genes)
    print('1111111')
    print(rootdir)
    if gene_num==1:
        file_path1,selects1 = single_gene_select_spacer(input_genes,Time)
        file_path2,selects2  = single_gene_mpcp(input_genes,Time)
    elif gene_num==2:
        file_path,select_spacers=two_gene_select_spacer(input_genes,Time)
    elif gene_num>2:
        file_path = multi_gene_select_spacer(input_genes,Time)


    #input_genes='ENSG00000002834'
    #select_spacer, select_spacer_ex, on_site_chr, on_site_strand, on_site_pos, on_site_cut_eff, on_site_diff, off_site_num, off_site_seq, mpcp, mpcp_chr, mpcp_strand, mpcp_pos, mpcp_cut_eff, mpcp_mis_num, mpcp_diff=single_gene_mpcp(gene)
    #single_gene_select_spacer(input_genes,111)
    #input_genes = 'ENSG00000002834,ENSG00000005073,ENSG00000005339'
    #input_genes = 'ENSG00000005073,ENSG00000002834'
    #file_path=two_gene_select_spacer(input_genes,131)
    #file_path=multi_gene_select_spacer(input_genes,121)
