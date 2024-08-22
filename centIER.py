#!/usr/bin/env python
import argparse
import sys,subprocess,os
import pyfastx,random
from collections import defaultdict
import numpy as np
import pandas as pd
import shutil
from scipy.sparse import coo_matrix
from translate_seq import six_frame_translate

script_path = os.path.dirname(os.path.abspath(__file__)) ### find the centrie path
Version = "V_3.0.1"



def argparseFunc():
    software = "  The software is designed to predict the location of centromeric \
regions in genomes assembled using T2T technology. It enables users to identify and annotate \
sequences within these regions.\nif you use this program, please cite the atricle:\n \t \
Xu, D. et al. CentIER: accurate centromere identification for plant genome. Plant Communications 101046 (2024) doi:10.1016/j.xplc.2024.101046."
    parser = argparse.ArgumentParser(description=software)
    parser.add_argument("genome",help="Genome fasta file need to annotation")
    parser.add_argument('-o',"--output", default='./CentIER_final_results', help='output path')
    parser.add_argument("--gff",help="optional, annotation file (gff or gtf format)")
    parser.add_argument("-k","--kmer_size",type=int,default=21,help="the size of kmer")
    parser.add_argument("-c","--center_tolerance",type=int,default=15,help="the fine-tuning size between two regions")
    parser.add_argument("--step_len",type=int,default=10000,help="the size between two regions")
    parser.add_argument('--mul_cents',action="store_true",help="This parameter is used to retain all potential centromeric regions, especially when the chromosomes of the species under study may not have a single, defined centromere; it should be included in such cases.")
    parser.add_argument('--matrix1', type=str, help='Path to matrix of 100000 file')
    parser.add_argument('--matrix2', type=str, help='Path to matrix of 200000 file')
    parser.add_argument('--bed1', type=str, help='Path to bed file of matrix1 ')
    parser.add_argument('--bed2', type=str, help='Path to bed file of matrix2')
    parser.add_argument('--MINGAP', type=int, default=2, help='Minimum gap value n*100000 (default: 2)')
    parser.add_argument('--SIGNAL_THRESHOLD', type=float, default=0.7, help='Signal threshold value (default: 0.7)')
    return parser.parse_args()

def get_interval(buck,name,threshold):
    centromeres = []
    dir_range={}
    temp = []
    for key in buck.keys():
        if(key < threshold):
            temp.extend(buck[key])
    temp.sort()
    interval = [0,0]
    intervals = []
    last_cordinate = 0
    for cordinate in temp:
        if(abs(cordinate - last_cordinate) <= step_len * center_tolerance):
            if(interval[1] == last_cordinate):
                if(last_cordinate == 0):
                    interval[0] = cordinate - step_len
                interval[1] = cordinate
            else:
                if(interval[1] - interval[0] > step_len):
                    intervals.append(list(interval))
                interval[0] = last_cordinate - step_len
                interval[1] = cordinate
        last_cordinate = cordinate
    if(interval[1] - interval[0] > step_len and interval not in intervals):
        intervals.append(list(interval))
    if(len(intervals) > 0):
        intervals.sort(key = lambda interval : interval[0])
        # print(intervals)
        for interval in intervals:
            centromeres.append([max(0,interval[0]),min(interval[1],chr_length[name])])
            # print([max(0,interval[0]),min(interval[1],len(fasta_sequence[name]))])
    dir_range[name]=centromeres
    return dir_range

def kmer_cal(file):
    # fasta_sequence={name:seq for name,seq in pyfastx.Fasta(file,build_index=False)}
    count=0
    for name,sequence in fasta_sequence.items():
        sequence=sequence.upper()
        chrid_list.append(name)
        total_length=len(sequence)
        chr_length[name]=total_length
        array = [0.0 for _ in range(step_len // 1000 + 1)]
        for i in range(len(sequence) - kmer_size + 1):
            temp_kmer = sequence[i : i + kmer_size]
            kmer[temp_kmer] += 1
            all_kmer[temp_kmer] += 1
            count += 1
            if(count % step_len == 0):
                length = len(kmer)
                array[length // 1000] += 1
                buck[length // 1000].append(count)
                kmer.clear()
        for i in range(len(array)):
            array[i] /= total_length // step_len + 1
        temp = 0
        threshold = len(array) - 1
        for i in range(len(array) - 1,0,-1):
            temp += array[i]
            if(temp >= 0.8):
                threshold = i
                break
        single=get_interval(buck,name,threshold)
        if type(single)==dict:arange.update(single)
        else:arange.update({name:''})
        kmer.clear()
        all_kmer.clear()
        buck.clear()
        count = 0
        for i in range(int(total_length/200000)+1):
            seq=sequence[i*200000:(i+1)*200000];std=i*200000;end=(i+1)*200000
            diseq={seq[i:i+20] for i in range(len(seq)-19)}#查找kmer的种类数
            if name not in kmer_count:kmer_count[name]=set()
            if len(diseq)<100000:
                kmer_count[name].add(std);kmer_count[name].add(end)
    return arange
    
def merge_regions(region):#单独处理每个得到的区间，便于有交集的区间合并
    final_result={};target_list=[]
    for i,all in region.items():
        length=len(all);all=sorted(all,key=lambda x:x[0])
        if i not in final_result:final_result[i]=[]
        for index,j in enumerate(all):
            if index<=length-2:
                if target_list==[]:target_list=j
                m1,n1=all[index+1]
                m,n=target_list
                if (n<m1 or m>n1)==False and index<length-2:target_list=[min(m,n,m1,n1),max(m,n,m1,n1)]
                elif index==length-2 or (n<m1 or m>n1)==True:
                    if (n<m1 or m>n1)==False:target_list=[min(m,n,m1,n1),max(m,n,m1,n1)]
                    final_result[i].append(target_list)
                    target_list=[]
            if index==length-1 and final_result[i]!=[]:
                m,n=final_result[i][-1]
                m1,n1=all[-1]
                if (n<m1 or m>n1)==False:
                    target_list=[min(m,n,m1,n1),max(m,n,m1,n1)]
                    if target_list not in final_result[i]:
                        final_result[i].append(target_list)
                else:
                    final_result[i].append(all[-1])
    return final_result

def ex_gff(file):
    with open(file) as f:
        count_dir={}
        gff_gene_number = {}
        for line in f:
            if '#' not in line and 'mRNA' in line:
                linelist=line.split()
                chrid=linelist[0];pos=int((int(linelist[3])+int(linelist[4]))*0.5/1000000)
                if chrid not in count_dir:count_dir[chrid]={}
                if pos not in count_dir[chrid]:count_dir[chrid][pos]=1
                else:count_dir[chrid][pos]+=1
        for j,count in count_dir.items():
            gff_gene_number[j]={}
            for i in range(max(count)):
                if i not in count:
                    gff_gene_number[j][str(i)+'-'+str(i+1)]=0
                else:
                    content=count[i]
                    gff_gene_number[j][str(i)+'-'+str(i+1)]=content
    return gff_gene_number

def creat_matrix(input):
    data = np.loadtxt(input, delimiter="\t", dtype=int)
    sparse_matrix = coo_matrix((data[:, 2], (data[:, 0], data[:, 1])))
    dense_matrix = sparse_matrix.toarray()
    return dense_matrix

def findGap(matrix):
    # This function is finding Candidates for Centremeres and they will be return
    # with dict.
    degrees = []
    # Calculate degrees of nodes
    for i in range(4, len(matrix)-3):
        ddq = 0
        for j in range(i-3, i):
            ddq += matrix[j][i]
        for j in range(i,i+3):
            ddq += matrix[i][j]
        degrees.append(ddq)
    avg = np.mean(sorted(degrees)[:-1000])

    variance = []
    for i in degrees:
        variance.append((i-avg)/avg)
    gap = {}
    errToler = {}
    singal = None
    for i in range(0, len(variance)):
        if singal and (-variance[i] > SIGNAL_THRESHOLD):
            gap [singal] += 1
        if (singal == None) and (-variance[i] > SIGNAL_THRESHOLD):
            gap [i] = 1
            singal = i
            errToler[i] = 0
        if singal and (-variance[i] <= SIGNAL_THRESHOLD):
            if errToler[singal] < 2:
                errToler[singal] += 1
                gap [singal] += 1
            else:
                k = i-1
                while -variance[k] <= SIGNAL_THRESHOLD:
                    gap[singal] -= 1
                    errToler[singal] -= 1
                    k -= 1
                singal = None
    candidates = {key+4: value for key, value in gap.items()}
    # Filter out the too short candidated gaps that are likely to be noise. 
    filter_candidates = {key: value for key, value in 
                         candidates.items() if value > MINGAP}
    return filter_candidates

def read_finer_matrix(input, bed):
    
    # get diag weight and return it as  list format
    finer_diag = np.zeros(len(bed.index)+1)
    with open(input, 'r') as file:
        i = 1
        for line in file:
            arr = line.strip().split("\t")
            if arr[0] == arr[1]:
                finer_diag[int(arr[0])] = arr[2]
    return finer_diag

def findPosition(gap, selfweights, bed1, bed2):
    centremere = {}
    mean_weight = np.mean(sorted(selfweights)[:-1000]) 
    for key, item in gap.items():
        if bed1[bed1["id"] == key]["chr"].iloc[0] != bed1[bed1["id"] == (key+item)]["chr"].iloc[0]:
            continue
        chr = bed1[bed1["id"] == key]["chr"].iloc[0]
        start = bed1[bed1["id"] == key]["start"].iloc[0]
        end = bed1[bed1["id"] == (key+item)]["start"].iloc[0]

        Lstart = int(bed2[(bed2["chr"] == chr) & (bed2["start"]==start)]["id"].iloc[0])
        Rstart = int(bed2[(bed2["chr"] == chr) & (bed2["start"]==end)]["id"].iloc[0])
        lbound, rbound = None, None
        for i in range(Lstart, Lstart + 10):
            if selfweights[i]/mean_weight < 0.4:
                lbound = bed2[bed2["id"]==i]["start"].iloc[0]
                break
        for i in range(Rstart, Rstart+10):
            if selfweights[i]/mean_weight > 0.7:
                rbound = bed2[bed2["id"]==i]["end"].iloc[0]
                break
        if chr not in centremere:
            centremere[chr] = []
        centremere[chr].append([lbound, rbound])
    return centremere

def strict_range(bat):
    with open(bat) as f:#对于单体长度，后续需要设置选项供用户选择
        monomer={}
        for line in f:
            if 'Sequence: ' in line:
                chrid=line.split()[1]
                monomer[chrid]={}
            if line.strip().count(' ')==14:
                linelist=line.split(' ')
                if 20<=int(linelist[2])<=1000 and float(linelist[3])>=10:
                    monomer[chrid][(int(linelist[0]),int(linelist[1]),linelist[13])]=float(linelist[3])
    monomer={i:sorted(j.items(),key = lambda x:x[1],reverse = True)[:10] for i,j in monomer.items()}#((14841146, 17128917, 'AGGCGTAAGAATTGTATCCTTGTTAAAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATCTCATGTGTATGATTGAGTATAAGAACTTAAACCGCAACCGGATCTTAA'), 12867.4)
    return monomer

def strict_range_repeat(dirfinal,bat):
    with open(bat) as f:#对于单体长度，后续需要设置选项供用户选择
        monomer={}
        for line in f:
            if 'Sequence: ' in line:
                chrid=line.split()[1]
                monomer[chrid]={}
            if line.strip().count(' ')==14:
                linelist=line.split(' ')
                if 100<=int(linelist[2])<=200 and float(linelist[3])>=10:
                    monomer[chrid][linelist[0]+' '+linelist[1]+' '+linelist[13]]=float(linelist[3])
    monomer={i:sorted(j.items(),key = lambda x:x[1],reverse = True) for i,j in monomer.items()}
    finalpredict={}
    for i,j in monomer.items():
        finalpredict[i]=''
        for m in j:
            infor=m[0].split()
            s=int(infor[0]);e=int(infor[1]);seq=infor[2]
            if i in dirfinal:
                for n in dirfinal[i]:
                    if n!=[]:
                        minvalue=n[0];maxvalue=n[-1]
                        if maxvalue<s or minvalue>e:
                            continue
                        else:
                            fmin=min([maxvalue,minvalue,s,e]);fmax=max([maxvalue,minvalue,s,e])
                            if fmax-fmin<500000:
                                fmax+=500000;fmin-=500000
                            elif 500000<=fmax-fmin<1000000:fmax+=300000;fmin-=300000
                            # finalpredict[i]=(min([maxvalue,minvalue,s,e]),max([maxvalue,minvalue,s,e]))
                            finalpredict[i]=(fmin,fmax,seq)
                            break
            if finalpredict[i]!='':break
    return finalpredict

def identify_period(content):
    period_copies={}
    for i in content:
        linelist=i.split()
        period_size=len(linelist[13]);copies=float(linelist[3])
        if period_size>=100 and copies>=2:
            if period_size not in period_copies:period_copies[period_size]=0
            period_copies[period_size]+=copies
    all_copy=list(period_copies.values())
    all_copy.sort(reverse=True)
    target_cop=set(all_copy[:5])
    target_period=[i for i,j in period_copies.items() if j in target_cop]
    return target_period

def identify_period_second(content,period):
    period_copies={}
    for i in content:
        linelist=i.split()
        period_size=len(linelist[13]);copies=float(linelist[3])
        r1=period-50;r2=period+50
        if r1<=0:r1=30
        if r2>=period_size>=r1 and copies>=2:
            if period_size not in period_copies:period_copies[period_size]=0
            period_copies[period_size]+=copies
    all_copy=list(period_copies.values())
    all_copy.sort(reverse=True)
    target_cop=all_copy[0]
    target_period=[i for i,j in period_copies.items() if j == target_cop]
    return target_period[0]

def search_range(content,f_period):
    copies_monomer={};monomer_range={};monomer_consensus={};size_range=[];cop_range={};total_range=[]
    consensus_max=max([int(line.split()[4]) for line in content if len(line.split()[13])==f_period])
    for line in content:
        linelist=line.split()
        st=int(linelist[0]);ed=int(linelist[1]);period_size=len(linelist[13]);copies=float(linelist[3]);consensus_size=int(linelist[4]);monomer=linelist[13]
        if period_size==f_period==consensus_size==consensus_max:
            size_range.append(copies)
            if copies not in cop_range:cop_range[copies]=[]
            cop_range[copies].append(st);cop_range[copies].append(ed)
            copies_monomer[copies]=monomer
    size_range.sort(reverse=True)
    top1=size_range[0]
    # for i in top5:
    for j in cop_range[top1]:total_range.append(j)
    monomerseq=copies_monomer[top1]
    minr=min(total_range);maxr=max(total_range)
    if maxr-minr>2000000:
        maxr-=300000;minr+=300000
    if maxr-minr<1000000:
        maxr+=300000;minr-=300000
        if minr<0:minr+=300000
    minr=min(total_range);maxr=max(total_range)
    return [monomerseq,minr,maxr]

def find_enrichment_with_bin_size(big_interval, small_intervals, bin_size=1):
    start_big, end_big = big_interval
    num_bins = ((end_big - start_big) // bin_size) + 1
    enrichment_map = defaultdict(int)
    for start_small, end_small in small_intervals:
        start_small = max(start_small, start_big)
        end_small = min(end_small, end_big)
        if start_small <= end_small:
            for bin_index in range((start_small - start_big) // bin_size, (end_small - start_big) // bin_size + 1):
                if bin_index < num_bins:
                    enrichment_map[bin_index] += 1
    return dict(enrichment_map)

def search_ltr1(seqfile):
    exe=os.path.abspath(script_path+'/bin/ltr_finder/ltr_finder')
    arg=[exe,'-D','20000','-d','1000','-L','3500','-l','100','-p','20','-C','-M','0.9',seqfile,'>',seqfile+'_l1.txt']
    result=subprocess.Popen(' '.join(arg),shell=True)
    result.wait()
    ltr1={}
    with open(seqfile+"_l1.txt") as f:
        for line in f:
            if ">Sequence:" in line:
                ID=line.split()[1].replace("_centromere_seq","")
                ltr1[ID]=""
            elif "Location : " in line:
                pos=line.split()[2]+"\t"+line.split()[4]
                ltr1[ID]=pos
    return ltr1

def search_ltr2(file):
    fname=os.path.basename(file)
    if os.path.exists('./gtltr'):
        shutil.rmtree('./gtltr')
        os.mkdir('./gtltr')
    else:
        os.mkdir('./gtltr')
    source=file;seqfile="./gtltr/"+fname
    shutil.copy(source, seqfile)
    exe='gt'
    arg=[exe,'suffixerator','-db',seqfile,'-indexname',seqfile,'-tis','-suf','-lcp','-des','-ssp','-sds','-dna']
    step1=subprocess.Popen(arg)
    step1.wait()
    arg2=[exe,'ltrharvest','-index',seqfile,'>','./gtltr/genometools_result.txt']
    step2=subprocess.Popen(' '.join(arg2),shell=True)
    step2.wait()
    with open("./gtltr/genometools_result.txt") as f:
        ltr2={}
        for line in f:
            if "#" not in line:
                linelist=line.split()
                s,e,ID=linelist[0],linelist[1],int(linelist[-1])
                if ID not in ltr2:ltr2[ID]=[]
                ltr2[ID].append((s,e))
    return ltr2

def exseq_translate(prefix,fa_file,script_path):
    fa=pyfastx.Fasta(fa_file)
    with open(prefix+"_ltr_position.txt") as f,open(prefix+"_ltr_position.txt_seq.txt","w") as f1:
        for line in f:
            linelist=line.split()
            chrid=linelist[0];ranges=(int(linelist[1]),int(linelist[2]))
            seq=fa.fetch(chrid,ranges)
            f1.write(">"+"*".join(linelist)+"\n"+seq+"\n")

    with open(prefix+"_ltr_position.txt_seqt.txt","w") as f:
        inSeq=prefix+"_ltr_position.txt_seq.txt"
        six_frame_translate(inSeq, f)
    filepath=prefix+"_ltr_position.txt_seqt.txt"
    hmmsearch= script_path + '/bin/hmmsearch' 
    args=[' ','--noali','--tblout',prefix+'LTR-hmmresult.txt','-E','1e-5',script_path+'/bin/REXdb.hmm',filepath]
    args=hmmsearch+" ".join(args)
    result=subprocess.Popen(args,shell=True) #这里，kegg的结果文件"./pfam-hmmresult.txt"应该存放在本目录下，用完后删除即可
    result.wait()
    with open(prefix+'LTR-hmmresult.txt') as f, open(prefix+"_LTR_positionb.txt","w") as f1:
        content={line.split()[0].split("|")[0].replace("*","\t")+"\t"+line.split()[2]+"\n" for line in f if "#" not in line}
        f1.write("".join(content))

def draw_result(prefix):
    with open(prefix+"_centromere_range.txt") as f:
        cen_length={line.split()[0]:int(line.split()[2])-int(line.split()[1]) for line in f}#着丝粒长度
    with open(prefix+"_monomer_in_centromere.txt") as f:#monomer在着丝粒上的位置
        repeat_position={line.split()[0]:line.split()[1]+"\t"+line.split()[2] for line in f}
    with open(prefix+"_LTR_positionb.txt") as f:
        ltr_monomer={};monomer_type=set();compare=set()
        for line in f:
            linelist=line.split()
            chrid=linelist[0];stp=int(linelist[1]);end=int(linelist[2]);length=str(abs(end-stp));monomer_ID=linelist[3]
            monomer_type.add(monomer_ID);compare.add(str(stp)+"\t"+length)
            if chrid not in ltr_monomer:ltr_monomer[chrid]=[]
            ltr_monomer[chrid].append(linelist[1]+"\t"+length+"\t"+monomer_ID)
    with open(prefix+"_ltr_position.txt") as f:
        rt_position={}
        for line in f:
            linelist=line.split()
            ID=linelist[0];p1=linelist[1];p2=linelist[2];p3=str(int(p2)-int(p1))
            if ID not in rt_position:rt_position[ID]=[]
            if p1+"\t"+p3 not in compare:rt_position[ID].append(p1+"\t"+p3)
    ql=max(cen_length.values())
    if 100<ql*10**-6<1000:unit=10**-6;unit_text="(Mb)"
    elif 100<ql*10**-5<1000:unit=10**-5;unit_text="(100 Kb)"
    elif 100<ql*10**-4<1000:unit=10**-4;unit_text="(10 Kb)"
    elif 100<ql*10**-3<1000:unit=10**-3;unit_text="(Kb)"
    elif 100<ql*10**-2<1000:unit=10**-2;unit_text="(100 bp)"
    elif 100<ql*10**-1<1000:unit=10**-1;unit_text="(10 bp)"
    elif 100<ql<1000:unit=1;unit_text="(bp)"
    chr_line='    <path d="M {},{} h {}"/>'
    kmer_mark='    <path d="M {},{} v {}"/>'
    chr_text='    <text x="{}" y="{}">{}</text>'
    x_range=ql*unit
    n=0;chr_seg={};r=lambda:random.randint(0,255);repc=[];ltr_c=[];mark_line=[];mark_text=[];chr_tests=[];unit_test=[];ltr_monomers={}
    monomer_type_color={i:"#%02X%02X%02X"%(r(),r(),r()) for i in monomer_type}
    for ID,length in cen_length.items():
        n+=1
        query_length=length*unit
        xp=5;yp=10*n
        chr_color="#%02X%02X%02X"%(r(),r(),r())
        query_chr=chr_line.format(xp,yp,query_length)
        chr_seg[query_chr]=chr_color
        chr_tests.append(chr_text.format(2.5,yp,"CEN-"+ID))
        if ID in repeat_position:
            rp1=5+int(repeat_position[ID].split("\t")[0])*unit;rp2=abs((int(repeat_position[ID].split("\t")[1])-int(repeat_position[ID].split("\t")[0]))*unit)
            repc.append(chr_line.format(rp1,yp-3,rp2))#需要加标题
        if ID in rt_position:
            rt_lists=rt_position[ID]
            for line in rt_lists:
                ltr_c.append(chr_line.format(int(line.split("\t")[0])*unit+5,yp+2,int(line.split("\t")[1])*unit))
        if ID in ltr_monomer:
            ltr_monomerl=ltr_monomer[ID]
            monomertype=set();monomernumber=[]
            for index,i in enumerate(ltr_monomerl):
                lis=i.split()
                stp=int(lis[0])*unit+5
                if stp not in monomertype:monomertype.add(stp)
                else:
                    monomernumber.append(stp)
                    stp=stp+0.02*monomernumber.count(stp)
                lth=int(lis[1])*unit;ID=lis[2]
                if ID not in ltr_monomers:ltr_monomers[ID]=[]
                ltr_monomers[ID].append(chr_line.format(stp,yp,lth))
        x1,y1,length=5,yp,query_length
        for k in range(int((length)/0.2)):
            x3=x1+k*0.2
            if k%10==0:
                lm=kmer_mark.format(x3,y1-1,0.5)#大刻度
                mark_line.append(lm)
                markt=chr_text.format(x3,y1-1.1,str(int(k/5)))#只有在大刻度上才有数字数字
                mark_text.append(markt)
            else:
                lm=kmer_mark.format(x3,y1-0.75,0.25)#小刻度
                mark_line.append(lm)
        unit_test.append(chr_text.format(x3+1.5,y1-1.1,unit_text))
    with open(prefix+"_draw_cen.svg","w") as f:
        x=x_range+10;y=len(cen_length)*10+10
        max_i=int((x-5)/6);l=len(monomer_type_color)
        xlength=len(monomer_type_color)*6+5
        if xlength>x:y=len(cen_length)*10+10*(int((xlength-x)/x)+2)
        if y<y-3.5+int((l+1)/max_i)*3:
            file_title='<svg width="{}px" height="{}px" xmlns="http://www.w3.org/2000/svg" version="1.1">\n'.format(x,y-3+int((l+1)/max_i)*3)
        else:
            file_title='<svg width="{}px" height="{}px" xmlns="http://www.w3.org/2000/svg" version="1.1">\n'.format(x,y)
        mark_l='  <g id="mark" stroke="#8e8d8d" stroke-width="0.02">\n'+"\n".join(mark_line)+"\n  </g>\n"
        markfont='  <g id="mark_text" fill="black" font-size="0.5" text-anchor="middle" font-family="Times New Roman">\n'+"\n".join(mark_text)+"\n  </g>\n"
        unitfont='  <g id="mark_text" fill="black" font-size="0.5" text-anchor="middle" font-family="Times New Roman">\n'+"\n".join(unit_test)+"\n  </g>\n"
        chr_ID='  <g id="CHR_text" fill="black" font-size="1" text-anchor="middle" font-family="Times New Roman">\n'+"\n".join(chr_tests)+"\n  </g>\n"
        query_segment=""
        for infor,color in chr_seg.items():
            single='  <g id="chrs" stroke="{}" stroke-width="1">\n'.format(color)+infor+"\n  </g>\n"
            query_segment+=single
        colorlista=["#8ecfc9","#ffbe7a","#fa7f6f","#82b0d2","#beb8dc","#2878b5","#9ac9db","#f8ac8c","#c82423","#ff8884","#14517c","#2f7fc1","#96c37d","#f3d266","#d8383a","#a9b8c6","#c497b2","#8e8bfe","#fe99a2","#934b43","#d76364","#ef7a6d","#63e398","#b1ce46","#f1d77e","#9394e7","#5f97d2","#9dc3e7","#a1a9d0","#f0988c","#b883d3","#c4a5de","#f6cae5","#96cccb"]
        ltr_color=random.choice(colorlista)
        # ltr_c.append(chr_line.format(5+len(monomer_type_color)*6,y-3,2))
        monomer_legend={ID:chr_line.format(3+6*(i%max_i),y-5+int(i/max_i)*3,2) for i,ID in enumerate(monomer_type_color)}
        monomer_legend_text={chr_text.format(4+6*(i%max_i),y-3.5+int(i/max_i)*3,ID) for i,ID in enumerate(monomer_type_color)}
        ltr_c.append(chr_line.format(3+6*(l%max_i),y-5+int(l/max_i)*3,2))
        monomer_legend_text.add(chr_text.format(4+6*(l%max_i),y-3.5+int(l/max_i)*3,"unknown LTR"))
        ltr_segments='  <g id="ltr" stroke="{}" stroke-width="0.5" opacity=".5">\n'.format(ltr_color)+"\n".join(ltr_c)+"\n  </g>\n"
        colorlista.remove(ltr_color)
        re_color=random.choice(colorlista)
        l=l+1
        monomer_legend_text.add(chr_text.format(4+6*(l%max_i),y-3.5+int(l/max_i)*3,"monomer"))
        repc.append(chr_line.format(3+6*(l%max_i),y-5+int(l/max_i)*3,2))
        re_segments='  <g id="re" stroke="{}" stroke-width="0.5">\n'.format(re_color)+"\n".join(repc)+"\n  </g>\n"
        monomer_legends=""
        monomer_ID='  <g id="CHR_text" fill="black" font-size="1" text-anchor="middle" font-family="Times New Roman">\n'+"\n".join(monomer_legend_text)+"\n  </g>\n"
        for ID,infor in monomer_legend.items():
            a='  <g id="monomer_legend" stroke="{}" stroke-width="1" stroke-opacity="0.5">\n'.format(monomer_type_color[ID])+infor+"\n  </g>\n"
            monomer_legends+=a
        ltr_monomer_results=""
        if ltr_monomers!={}:
            for ID,infor in ltr_monomers.items():
                ltr_monomer_color=monomer_type_color[ID]
                ltr_monomer_result='  <g id="chrs" stroke="{}" stroke-width="1.2" stroke-opacity="0.5">\n'.format(ltr_monomer_color)+"\n".join(infor)+"\n  </g>\n"
                ltr_monomer_results+=ltr_monomer_result
        end='  </svg>'
        f.write(file_title+mark_l+markfont+query_segment+ltr_segments+re_segments+ltr_monomer_results+chr_ID+unitfont+monomer_legends+monomer_ID+end)

def statistics_ltr(prefix):
    with open(prefix+"_LTR_positionb.txt") as f,open(prefix+"_LTR_statistics.txt","w") as f1:
        dirs={};alltype=[]
        for line in f:
            linelist=line.split()
            chrid=linelist[0];types=linelist[-1]
            if chrid not in dirs:
                dirs[chrid]=[]
            dirs[chrid].append(types)
        for i in dirs.values():
            alltype+=i
        alltype=list(set(alltype))
        alltype.sort()
        alltype.insert(0,"ID")
        f1.write("\t".join(alltype)+"\n")
        for ID,content in dirs.items():
            infor=[str(content.count(i)) for index,i in enumerate(alltype) if index>=1]
            infor.insert(0,ID)
            f1.write("\t".join(infor)+"\n")

if __name__ == '__main__':
    trf, gt, ltr = 0,0,0
    RED = '\033[91m'
    RESET = '\033[0m' ### give a error color
    args = argparseFunc()
    fasta = args.genome
    gff=args.gff
    print(gff)
    kmer_size = args.kmer_size
    center_tolerance = args.center_tolerance
    step_len = args.step_len
    MINGAP = args.MINGAP
    SIGNAL_THRESHOLD = args.SIGNAL_THRESHOLD
    kmer = defaultdict(int)
    all_kmer = defaultdict(int)
    buck = defaultdict(list)
    if not os.path.exists(args.output):
        os.mkdir(args.output)
    else:
        print(f"Folder '{args.output}' already exists. Skipping creation.")
    output = args.output+"/"

    if args.mul_cents:
        mul_cents=True
    else:
        mul_cents=False
    kmer_count={};listl=[];ranges={};dirfinal={} #chr_length_list记录的是染色体ID和长度
    chrid_list=[];chr_length={};arange={};fasta_sequence={};repeat_file=False
    fasta_sequence={name:seq for name,seq in pyfastx.Fasta(fasta,build_index=False)}
    #The Linux system should be a 64-bit system.
    # try:
    print ('Searching for tandem repeat sequences')
    prefix = os.path.basename(fasta) ### Temporary file prefix


    if trf:
        arg=[script_path+"/bin/trf409.linux64",fasta,'2','5','7','80','10','50','2000','-h','-f','-d','-m','-l','15']
        TRF_search=subprocess.Popen(arg)
        print ('searching for all LTRs of the input genome')

    database = prefix+'_LTR_database'
    
    if gt:
     ### database put in tmp path
        a=shutil.which('gt')
        if a!=None:
            exe=a
            arg=[exe,'suffixerator','-db',fasta,"-indexname", database,'-tis','-suf','-lcp','-des','-ssp','-sds','-dna']
            step1=subprocess.Popen(' '.join(arg),shell=True)
            step1.wait()
            arg=[exe,'ltrharvest','-index',database,'-similar','90','-vic','10','-seed','20','-seqids','yes','-minlenltr','100','-maxlenltr','7000','-mintsd','4','-maxtsd','6','-motif','TGCA','-motifmis','1','>',database+'.harvest.scn']
            step2=subprocess.Popen(' '.join(arg),shell=True)
            exe=script_path+'/bin/ltr_finder/ltr_finder'
            arg=[exe,'-D','15000','-d','1000','-L','7000','-l','100','-p','20','-C','-M','0.9',fasta,'>', database+'.finder.scn']
            step3=subprocess.Popen(' '.join(arg),shell=True)
        else:
            sys.stderr.write(f"{RED}The gt software has not been detected in the \
                            environment variables of the server. Please download \
                            and install it, and then add it to the environment \
                            variables. The download link is \
                            https://github.com/genometools/genometools..{RESET}\n")
            sys.exit(1)
    
    arange=kmer_cal(fasta)
    # print('fasta_seq',fasta_sequence)
    kmer_dir=merge_regions(arange)
    file_size=os.path.getsize(fasta)
    if file_size <=500*1000000:
        repeat_file=True
        for chr,infor in kmer_count.items():
            ranges[chr]=[]
            infor=list(infor)
            infor.sort()
            if infor!=[] and len(infor)>1:
                for index,content in enumerate(infor):
                    if content>=1000000:
                        if index==0:
                            listl.append(content)
                        else:
                            if infor[index]-infor[index-1]<=400000:
                                listl.append(content)
                            else:
                                ranges[chr].append(listl)
                                listl=[]
                                listl.append(content)
                listl=[]
        dirfinal={}
        for i,j in ranges.items():
            dirfinal[i]=[]
            for m in j:
                if len(m)>1 and max(m)-min(m)<5000000:
                    dirfinal[i].append(m)
    if gff:gff_gene_number=ex_gff(gff)
    if args.bed1 and args.bed2 and args.matrix1 and args.matrix2:
        bed1 = pd.read_csv(args.bed1, sep = "\t", header= None, 
                        names=['chr', 'start', 'end', 'id'])
        bed2 = pd.read_csv(args.bed2, sep = "\t", header= None,
                        names=['chr', 'start', 'end', 'id'])
        gap = findGap(creat_matrix(args.matrix1))
        finer_diag = read_finer_matrix(args.matrix2, bed2)
        hic_dir = findPosition(gap, finer_diag, bed1, bed2) #hic_dir是收集hic信号的
        hic_dir=merge_regions(hic_dir)
    else:hic_dir={}
    
    if gt:
        step2.wait()
        step3.wait()
    
    
    if ltr:
        exe='LTR_retriever'
        a=shutil.which('LTR_retriever')
        if a!=None:
            # arg=exe+' -genome '+fasta+' -inharvest '+database+'.harvest.scn '+'-infinder '+database+'.finder.scn '+'-threads 60 -u 4.02e-9'
            arg=[exe,'-genome',fasta, '-inharvest',database+'.harvest.scn','-infinder',database+'.finder.scn','-threads','60','-u','4.02e-9']
            step4=subprocess.Popen(arg)
            step4.wait()
        else:
            sys.stderr.write(f"{RED}The LTR_retriever software has not been detected \
                            in the environment variables of the server. Please \
                            download and install it, and then add it to the environment \
                            variables. The download link is \
                            https://github.com/oushujun/LTR_retriever.{RESET}\n")
            sys.exit(1)
    if trf:
        TRF_search.wait()
    LTR_region={};LTR_level={};level_dir={}
    if os.path.exists(prefix+'.out.LTR.distribution.txt') or os.path.exists(prefix+'.mod.out.LTR.distribution.txt'):
        if os.path.exists(prefix+'.out.LTR.distribution.txt'):file=prefix+'.out.LTR.distribution.txt'
        else:file=prefix+'.mod.out.LTR.distribution.txt'
        with open(file) as f:
            content=f.readlines()[2:]
            if 'mod' in file:t=6
            else:t=4
            for line in content:
                linelist=line.split()
                chrid=linelist[0];level=float(linelist[t])
                if chrid not in LTR_level:LTR_level[chrid]=[]
                LTR_level[chrid].append(level)
            level_dir={i:set(sorted(j,reverse=True)[:10]) for i,j in LTR_level.items()}#原定[:3]
            for line in content:
                linelist=line.split()
                chrid=linelist[0];st=int(linelist[1]);ed=int(linelist[2]);level=float(linelist[t])
                if level in level_dir[chrid]:
                    if chrid not in LTR_region:LTR_region[chrid]=[]
                    LTR_region[chrid].append([st,ed])
        tf_dir=merge_regions(LTR_region)
    else:
        sys.stderr.write(f"{RED}Error: 'It appears that some errors occurred \
                         during the analysis of LTRs.You can contact the author \
                         to assist in resolving the issue.'{RESET}\n")
        sys.exit(1)

    bat=prefix+".2.5.7.80.10.50.2000.dat"
    if repeat_file==False or mul_cents==True:
        t_cont=[];k_cont=[];tf_cont=[];hic_cont=[];final_result={};target_list=[];number=0 #number用于存储交集的次数
        region_percent={};tandem_repeat={};number_list={};precise_range={}
        monomer_select=strict_range(bat)
        tandem_dir={}
        for i,j in monomer_select.items():
            tandem_dir[i]=[[m[0][0],m[0][1]] for m in j]
        tandem_dir2=merge_regions(tandem_dir)
        all_dir={}
        for i in chrid_list:
            all=[]
            # print(i)
            if i in tandem_dir:
                t_cont=tandem_dir[i];t_cont2=tandem_dir2[i]
                if len(t_cont2)/len(t_cont)<0.5:index=2
                else:index=1
                all+=index*t_cont
            if i in kmer_dir:k_cont=kmer_dir[i];all+=k_cont #;print('k_cont',k_cont)
            if i in tf_dir:tf_cont=tf_dir[i];all+=tf_cont #;print('tf_cont',tf_cont)
            if i not in final_result:final_result[i]=[]
            if i not in region_percent:region_percent[i]=[]
            if i not in number_list:number_list[i]=[]
            # print(i);print('t_cont',t_cont);print('k_cont',k_cont)
            if all!=[]:
                # print('all',all)
                all_dir[i]=all
                length=len(all);all=sorted(all,key=lambda x:x[0])
                # if i=='GWHBCKZ00000004':print('all',all)
                for index,j in enumerate(all):
                    if index<=length-2:
                        if target_list==[]:target_list=j
                        m1,n1=all[index+1]
                        m,n=target_list
                        if (n<m1 or m>n1)==False and index<length-2:target_list=[min(m,n,m1,n1),max(m,n,m1,n1)];number+=1
                        elif index==length-2 or (n<m1 or m>n1)==True:
                            if (n<m1 or m>n1)==False:target_list=[min(m,n,m1,n1),max(m,n,m1,n1)];number+=1
                            # final_result[i].append(target_list)
                            region_percent[i].append(list(target_list)+[number])
                            number_list[i].append(number)
                            target_list=[];number=0
                    if index==length-1:
                        m,n,l=region_percent[i][-1]
                        m1,n1=all[-1]
                        if (n<m1 or m>n1)==False:
                            target_list=[min(m,n,m1,n1),max(m,n,m1,n1)];number+=1
                            if target_list not in region_percent[i]:
                                # final_result[i].append(target_list)
                                region_percent[i].append(list(target_list)+[number])
                                number_list[i].append(number)
                        else:
                            # final_result[i].append(all[-1])
                            region_percent[i].append(list(all[-1])+[number])#第三项是可信度
                            number_list[i].append(number)
                        target_list=[]
        for i,content in region_percent.items():
            max_number=max(number_list[i])
            content_sort=sorted(content,key=lambda x:x[2],reverse=True)
            # print(content_sort)
            for j in content_sort:
                a,b,c=j
                if mul_cents==False and c==max_number:
                    final_result[i].append(j)
                elif mul_cents==True and c>=1: #预测多着丝粒情况
                    final_result[i].append(j)
        if mul_cents==False:
            s=0;e=0;all_list=[]
            bin_size = 500000
            lists=[];list_max=[]
            for i,content in final_result.items():
                chrlength=chr_length[i]
                all=all_dir[i]
                if len(content)>=2:
                    for r in content:
                        a,b,c=r
                        tf_cont=tf_dir[i]
                        for j in tf_cont:
                            d,e=j
                            if a<=e and b>=d:
                                big_interval=[a,b]
                                enrichment_regions = find_enrichment_with_bin_size(big_interval, all, bin_size)
                                confidence=sorted({j for i,j in enrichment_regions.items()},reverse=True)
                                thre=confidence[1];max_thre=confidence[0]
                                for bin_index, count in enrichment_regions.items():
                                    if count>=thre:
                                        start_of_bin = bin_index * bin_size + big_interval[0]
                                        end_of_bin = (bin_index + 1) * bin_size + big_interval[0] - 1
                                        if count==max_thre:list_max.append(start_of_bin);list_max.append(end_of_bin)
                                        lists.append(start_of_bin);lists.append(end_of_bin)
                                list_max=sorted(list_max)
                                max_v1=list_max[0];max_v2=list_max[-1]
                                s=min(lists);e=max(lists)
                                if s==0:
                                    if max_v1<e:s=max_v1
                                    else:s=1
                                if e-s>5000000:
                                    if e-max_v2>max_v1-s:e=max_v2
                                    else:
                                        s=max_v1
                                all_list.append((i,s,e))
                                lists=[]
                                list_max=[]
                                break
                            else:continue
                elif len(content)==1:
                    for r in content:
                        a,b,c=r
                        big_interval=[a,b]
                        enrichment_regions = find_enrichment_with_bin_size(big_interval, all, bin_size)
                        confidence=sorted({j for i,j in enrichment_regions.items()},reverse=True)
                        thre=confidence[1];max_thre=confidence[0]
                        if thre<2:thre==2
                        for bin_index, count in enrichment_regions.items():
                            if count>=thre:
                                start_of_bin = bin_index * bin_size + big_interval[0]
                                end_of_bin = (bin_index + 1) * bin_size + big_interval[0] - 1
                                if count==max_thre:list_max.append(start_of_bin);list_max.append(end_of_bin)
                                lists.append(start_of_bin);lists.append(end_of_bin)
                        list_max=sorted(list_max)
                        max_v1=list_max[0];max_v2=list_max[-1]
                        s=min(lists);e=max(lists)
                        if s==0:
                            if max_v1<e:s=max_v1
                            else:s=1
                        if e-s>5000000:
                            unit=int(((int((e-s)/1000000))*1000000)/6)
                            e-=unit;s+=unit
                        if e-s<1000000:s-=500000;e+=500000
                        if s<0:s=1
                        if e>chrlength:e=chrlength
                        all_list.append((i,s,e))
                        lists=[]
                        list_max=[]
            print('all_list',all_list)
            save=open(output + prefix+"_monomer_seq.txt","w")
            for j in all_list:
                i,s,e=j;detail=[]
                for d in monomer_select[i]:
                    a,b,c=d[0]
                    if a<s and b<e:detail.append((b-s,c))
                    elif s<=a and b<=e:detail.append((b-a,c))
                    elif s<a<e and b>e:detail.append((e-a,c))
                    elif a<=s and b>=e:detail.append((e-s,c))
                if detail!=[]:
                    detail_s=sorted(detail,key=lambda x:x[0],reverse=True)
                    save.write(i+'\t'+detail_s[0][1]+'\n')
            save.close()

        else:
            all_list=[]
            save=open(output + prefix +"_monomer_seq.txt","w")
            for i,j in final_result.items():
                for m in j:
                    s,e,c1=m;detail=[]
                    all_list.append((i,s,e))
                    for d in monomer_select[i]:
                        a,b,c=d[0]
                        if a<s and b<e:detail.append((b-s,c))
                        elif s<=a and b<=e:detail.append((b-a,c))
                        elif s<a<e and b>e:detail.append((e-a,c))
                        elif a<=s and b>=e:detail.append((e-s,c))
                    if detail!=[]:
                        detail_s=sorted(detail,key=lambda x:x[0],reverse=True)
                        save.write(i+'\t'+detail_s[0][1]+'\n')
            save.close()

    elif repeat_file==True and mul_cents==False:
        precise_range={};all_list=[]
        st_range=strict_range_repeat(dirfinal,bat)
        with open(fasta+".2.5.7.80.10.50.2000.dat") as f:
            allfor={};final_period=[];ID_period={};all_range=[]
            for line in f:
                if "Sequence: " in line:
                    ID=line.split()[1]
                    allfor[ID]=[]
                else:
                    if line.count(" ") == 14:allfor[ID].append(line)
            for ID,lists in allfor.items():
                p=identify_period(lists)
                all_range += p
            f_period=max(all_range,key = all_range.count)
            for ID,lists in allfor.items():
                p=identify_period_second(lists,f_period)
                ID_period[ID]=p
        print('precise_range1',precise_range)

        save=open(output + prefix+"_monomer_seq.txt","w")
        for ID,infors in allfor.items():
            strict_r=st_range[ID]
            if strict_r=='':
                seq,st,ed=search_range(infors,ID_period[ID])
                precise_range[ID]=(st,ed)
            else:
                seq=strict_r[2]
            save.write(ID+"\t"+seq+"\n")
        save.close()
        print('precise_range2',precise_range)
        repeat_candidate_range={}
        for i,j in arange.items():
            chrlength=chr_length[i]
            print('precise_range3',precise_range)
            b=0;range_list=[]
            if i in precise_range:s,e=precise_range[i]
            if i in st_range:strict_r=st_range[i]
            if strict_r!='':
                repeat_candidate_range[i]=(int(strict_r[0]),int(strict_r[1]))
                # rangefile.write(i+"\t"+str(strict_r[0])+"\t"+str(strict_r[1])+"\n")
            else:
                if e-s>2000000:
                    repeat_candidate_range[i]=(s,e)
                else:
                    if j!='':
                        jj={e2-s2:(s2,e2) for (s2,e2) in j if (s2>e or e2<s)==False}
                        if len(jj)!=0:
                            if len(jj)==1:
                                for k,v in jj.items():
                                    (s2,e2)=v
                            if len(jj)>1:
                                (s2,e2)=jj[max(jj)]
                            if s2<=s<e2<=e:
                                if e-s2>=3000000:
                                    s-=100000
                                    if s<=0:s=1
                                else:s=s2
                            elif s2<=s<e<=e2:
                                if e2-s2>=3000000:
                                    s-=100000
                                    if s<=0:s=1
                                else:s=s2;e=e2
                            elif s<=s2<e<=e2:
                                if e2-s>=3000000:
                                    e+=100000
                                    if e>chrlength:e=chrlength
                                else:e=e2
                            else:
                                s-=100000
                                if s<=0:s=1
                        else:
                            s-=100000
                            if s<=0:s=1
                        if e-s<500000:
                            e+=300000;s-=300000
                            if e>chrlength:e=chrlength
                            if s<=0:s=1
                        elif 500000<=e-s<1000000:
                            e+=200000;s-=200000
                            if e>chrlength:e=chrlength
                            if s<=0:s=1
                        repeat_candidate_range[i]=(s,e)
                    # rangefile.write(i+"\t"+str(s)+"\t"+str(e)+"\n")
        LTR_regiont={};LTR_number={};number=0
        for i,content in LTR_region.items():
            length=len(content);LTR_number[i]={}
            target_list=[]
            for index,j in enumerate(content):
                if index<=length-2:
                    if target_list==[]:target_list=j
                    m1,n1=content[index+1]
                    m,n=target_list
                    if (n<m1 or m>n1)==False and index<length-2:target_list=[min(m,n,m1,n1),max(m,n,m1,n1)];number+=1
                    elif index==length-2 or (n<m1 or m>n1)==True:
                        if (n<m1 or m>n1)==False:target_list=[min(m,n,m1,n1),max(m,n,m1,n1)];number+=1
                        if i not in LTR_regiont:LTR_regiont[i]=[]
                        LTR_regiont[i].append(target_list)
                        LTR_number[i][number]=target_list
                        target_list=[];number=0
        for ID,content in LTR_regiont.items():
            if ID in repeat_candidate_range:cs,ce=repeat_candidate_range[ID]
            else:cs=0;ce=0
            csl,cel=LTR_region[ID][-1]
            grn1=int(cs/1000000);grn2=int(ce/1000000)
            a=0
            LTRs_number=max(LTR_number[ID])
            for crange in content:
                s1,e1=crange
                if crange==LTR_number[ID][LTRs_number]:target_ltr=crange
                if s1<=cs<ce<=e1 or cs<=s1<e1<=ce:
                    lists=sorted([s1,e1,cs,ce])
                    all_list.append((ID,lists[1],lists[2]));a=1;break
                    # rangefile.write(ID+"\t"+str(lists[1])+"\t"+str(lists[2])+"\n");a=1;break
                else:
                    if cs<s1<ce:
                        all_list.append((ID,cs,ce));a=1;break
                    elif cs==ce==0:
                        all_list.append((ID,s1,e1));a=1;break
                        # rangefile.write(ID+"\t"+str(s1)+"\t"+str(e1)+"\n");a=1;break
                    elif s1>ce and s1-ce<=100000:
                        all_list.append((ID,cs,e1));a=1;break
                        # rangefile.write(ID+"\t"+str(cs)+"\t"+str(e1)+"\n");a=1;break
                    elif cs<e1<ce:
                        all_list.append((ID,s1,e1));a=1;break
                        # rangefile.write(ID+"\t"+str(s1)+"\t"+str(e1)+"\n");a=1;break
                    elif cs>e1 and cs-e1<=100000:
                        all_list.append((ID,s1,ce));a=1;break
                        # rangefile.write(ID+"\t"+str(s1)+"\t"+str(ce)+"\n");a=1;break
            if a==0:
                rep_gene_number=0;ltr_gene_number=0
                target=[]
                for tl in content:target+=tl
                min_value=min(target);max_value=max(target)
                if gff:
                    content=gff_gene_number[ID]
                    for i,j in content.items():
                        n1=int(i.split('-')[0])
                        if int(min_value/1000000)<=n1<=int(max_value/1000000):ltr_gene_number+=1
                        if grn1<=n1<=grn2:rep_gene_number+=1
                    if rep_gene_number<=ltr_gene_number:
                        all_list.append((ID,cs,ce))
                        # rangefile.write(ID+"\t"+str(cs)+"\t"+str(ce)+"\n")
                    else:
                        if max_value-min_value>3000000:
                            if ce<=min_value or cs>=max_value:
                                all_list.append((ID,min_value,max_value))
                                # rangefile.write(ID+"\t"+str(min_value)+"\t"+str(max_value)+"\n")
                            else:
                                all_list.append((ID,cs,ce))
                                # rangefile.write(ID+"\t"+str(cs)+"\t"+str(ce)+"\n")
                        else:
                            all_list.append((ID,min_value,max_value))
                            # rangefile.write(ID+"\t"+str(min_value)+"\t"+str(max_value)+"\n")
                else:
                    seqs=fasta_sequence[ID]
                    if cs<1000000 or chr_length[ID]-ce<1000000:
                        if len(content)==1:
                            all_list.append((ID,min_value,max_value))
                            # rangefile.write(ID+"\t"+str(min_value)+"\t"+str(max_value)+"\n")
                        else:
                            dirs={}
                            for tl in content:
                                m,n=tl
                                seq1=seqs[m:n]
                                diseq={seq1[i:i+20] for i in range(len(seq1)-19)}
                                dirs[len(diseq)]=tl
                            for i,content in dirs.items():
                                if i==min(dirs):
                                    all_list.append((ID,content[0],content[1]))
                                    # rangefile.write(ID+"\t"+str(content[0])+"\t"+str(content[1])+"\n")
                    else:
                        cs,ce=target_ltr
                        if cs<2000000:
                            all_list.append((ID,csl,cel))
                            # rangefile.write(ID+"\t"+str(csl)+"\t"+str(cel)+"\n")
                        else:
                            all_list.append((ID,cs,ce))
    #     monomer_position=open(fasta+"_monomer_in_centromere.txt","w")
    #     with open(fasta+"_centromere_range.txt") as f:
    #         for line in f:
    #             linelist=line.split()
    #             ID=linelist[0];s=int(linelist[1]);e=int(linelist[2])
    #             if ID in precise_range:s1,e1=precise_range[ID]
    #             else:s1,e1,seq1=st_range[ID]
    #             s2=s1-s;e2=e1-s
    #             if s2<0:s2=1
    #             if e2-s2>e-s:e2=e
    #             monomer_position.write(ID+"\t"+str(s2)+"\t"+str(e2)+"\n")
    #     monomer_position.close()
    # print('all_list2',all_list)
    print("****************************** Start defining the centromere range... *************************")
    rangefile=open(output + prefix+"_centromere_range.txt","w")
    [rangefile.write(i[0]+'\t'+str(i[1])+'\t'+str(i[2])+'\n') for i in all_list]
    rangefile.close()
    print("****************************** Centromere range determination completed *************************")
    all_cen_seq=open(output + prefix+"_all_centromere_seq.txt","w")
    fa=pyfastx.Fasta(fasta)
    with open(output + prefix+"_centromere_range.txt") as f:
        for line in f:
            linelist=line.split()
            ID=linelist[0]
            start=int(linelist[1]);end=int(linelist[2])
            # print(chr_length)
            chrlength=chr_length[ID]
            if end>chrlength:end=chrlength
            if start==0:start=1
            seq_range=(start,end)
            all_cen_seq.write(">"+ID+"_centromere_seq\n")
            all_cen_seq.write(fa.fetch(ID, seq_range)+"\n")
    all_cen_seq.close()
    # subprocess.Popen('cp '+fasta+"_all_centromere_seq.txt ./CentIER_final_results",shell=True)
    # if os.path.exists("./inter_doc")==False:
    #     creat_folder=subprocess.Popen("mkdir ./inter_doc",shell=True)
    #     creat_folder.wait()

    print("****************************** Annotating the centromere *************************")
    monomer_position=open(output + prefix +"_monomer_in_centromere.txt","w")
    with open(output + prefix+ "_centromere_range.txt") as f:
        for line in f:
            linelist=line.split()
            ID=linelist[0];s=int(linelist[1]);e=int(linelist[2])
            if ID in precise_range:s1,e1=precise_range[ID]
            else:s1,e1,seq1=st_range[ID]
            s2=s1-s;e2=e1-s
            if s2<0:s2=1
            if e2-s2>e-s:e2=e
            monomer_position.write(ID+"\t"+str(s2)+"\t"+str(e2)+"\n")
    monomer_position.close()

    # subprocess.Popen('cp '+fasta+"_monomer_in_centromere.txt ./CentIER_final_results",shell=True)
    ltr_result=open(output + prefix+"_ltr_position.txt", "w")
    ltr1=search_ltr1(output + prefix+"_all_centromere_seq.txt")
    [ltr_result.write(i+"\t"+j+"\n") for i,j in ltr1.items() if j!=""]

    ltr2=search_ltr2(output + prefix+"_all_centromere_seq.txt")
    for index, ID in enumerate(ltr1):
        if index in ltr2:
            for i in ltr2[index]:#如果中间有个染色体没有ltr这里可能需要处理
                s,e=i
                ltr_result.write(ID+"\t"+str(s)+"\t"+str(e)+"\n")
    ltr_result.close()


    exseq_translate(output + prefix, fasta, script_path)
    draw_result(output + prefix)
    statistics_ltr(output + prefix)
    print("*****************CentIER running finished!**************************")
    # subprocess.Popen('rm -rf ./gtltr')
    # subprocess.Popen('rm '+fasta+"_monomer_seq.txt")
    # subprocess.Popen('rm '+fasta+"_centromere_range.txt")
    # subprocess.Popen('rm '+fasta+"_all_centromere_seq.txt")
    # subprocess.Popen('rm '+fasta+"_monomer_in_centromere.txt")
    # if args.clear:
    #     subprocess.Popen('rm '+fasta+'_*')
    #     subprocess.Popen('rm '+fasta+'.*')


