
import os 
import sys
import zipfile

fastq1=list(map(str,sys.argv[1].strip().split(','))) 
fastq2=list(map(str,sys.argv[2].strip().split(',')))

project_name= sys.argv[3]

wdir=''
genome_ref_dir = '' 
rsem_ref_dir=''
raw_dir=''

#fastaq 파일이 match되지 않을 때 - > 시스템종료 
if (len(fastq1) != len(fastq2)):
        os.system('echo fastq file is not matched')
        sys.exit()
        
for i in range(len(fastq1)):
        fastq1[i] = fastq1[i][:-3]
        fastq2[i] = fastq2[i][:-3]
        
#make directory
#trimming 전후로 qc를 체크
os.chdir(wdir)
os.makedirs('./'+project_name,exist_ok=True)
os.chdir('./'+project_name)
path_list=['./qc','./qc2', './trimming', './mapping_sorting', 'counting']
for path in path_list:
        os.makedirs(path, exist_ok=True)

cdir = os.getcwd()

#qc
for i in range(len(fastq1)):
        os.system('gzip -d '+raw_dir+fastq1[i]+'.gz '+raw_dir+fastq2[i]+'.gz')
        os.system('fastqc -o '+cdir+'/qc/ -f fastq '+raw_dir+fastq1[i]+' '+raw_dir+fastq2[i])
        
os.chdir('./qc')#easy to move os.chdir command 
cdir = os.getcwd()

#fastqc 결과가 모여있는 폴더에서 _fastqc.zip파일을 열고 그 안에서 summary.txt에 해당하는 파일을 연다. 그 후 안에 있는 내용들을 모아 하나의 파일로 정리
files = [file for file in os.listdir(cdir) if file.endswith('_fastqc.zip')]
all_summary = []

for file in files:
        archive = zipfile.ZipFile(file,'r')
        members = archive.namelist()
        
        fname = [member for member in members if 'summary.txt' in member][0]
        data=archive.open(fname)
        
        for line in data:
                all_summary.append(line)
        data.close()
        archive.close()
        
with open('all_summary.txt', 'wb') as f: #summary.txt 가 바이너리 파일이라서 쓸 때도 이렇게 써야 한다 
        for content in all_summary:
                f.write(content)
        os.system('echo "./qc/all_summary.txt -> qc results of all fastqc file"')


#trimming        
for i in range(len(fastq1)):
        project = project_name.replace('#',str(i+1)) # '#'을 통해서 1번 2번 3번...환자를 구분한다 
        os.system('cutadapt -a AGATCGGAAGAGC -g AGATCGGAAGAGC -q 30 -m 20 -o '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' -p '+wdir+project_name+'/trimming/trimmed_'+fastq2[i]+' '+raw_dir+fastq1[i]+' '+raw_dir+fastq2[i])
        os.system('fastqc -o '+wdir+project_name+'/qc2/ -f fastq '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' '+wdir+project_name+'/trimming/trimmed_'+fastq2[i])
        os.system('gzip '+raw_dir+fastq1[i]+' '+raw_dir+fastq2[i])

#trimming 후 qc 진행 
os.chdir(wdir+project_name+'/qc2/')# '/'게 있어야 이동한대 
cdir=os.getcwd()

files = [file for file in os.listdir(cdir) if file.endswith('_fastqc.zip')]
all_summary = []

for file in files:
        archive = zipfile.ZipFile(file,'r')
        members = archive.namelist()
        fname = [member for member in members if 'summary.txt' in member][0]
        data=archive.open(fname)

        for line in data:
                all_summary.append(line)
        data.close()
        archive.close()
with open('after_trimming_all_summary.txt','wb') as f:
        for content in all_summary:
                f.write(content)
        os.system('echo "./qc2/after_trimming_all_summary.txt -> qc results of all trimmed_fasta_qc"')

#star로 mapping rsem으로 counting 
for i in range(len(fastq1)):
        project= project_name.replace('#',str(i+1))
        #mapping
        os.system('STAR --runMode alignReads --runThreadN 16 --outFilterMultimapNmax 10 --alignIntronMin 61 --alignIntronMax 265006 --genomeDir '+genome_ref_dir+' --readFilesIn '+wdir+project_name+'/trimming/trimmed_'+fastq1[i]+' '+wdir+project_name+'/trimming/trimmed_'+fastq2[i]+' --outFileNamePrefix '+wdir+project_name+'/mapping_sorting/star_mapsort_'+project+'_'+' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM')
        #counting
        os.system('rsem-calculate-expression -p 8 --alignments --paired-end --strandedness reverse --no-bam-output '+wdir+project_name+'/mapping_sorting/star_mapsort_'+project+'_Aligned.toTranscriptome.out.bam '+rsem_ref_dir+' '+wdir+project_name+'/counting/count_'+project)
        #최종 파일 이름 r에서 사용하기 쉽게 변경 
        os.chdir(wdir+project_name+'/counting/')
        result_name = 'count_'+project+'.genes.results'
        os.system('mv '+result_name+' '+result_name+'.txt')
