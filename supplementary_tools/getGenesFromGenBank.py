#!/usr/bin/python3
import sys,gzip,os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
#This script uses BioPython version 1.83

#Example on how to run the script:
#python3 getGenesFromGenBank.py gene_list.txt file1.gb file2.gb ...


#function merge_genbank (from GitHub links: https://github.com/kblin/merge-gbk-records)
def merge(records, length=10, spacer='n'):
    """Merge multiple SeqRecords into one, using a defined spacer

    :param records: Iterable containing SeqRecords to be merged
    :param length: Length of the spacer in kbp
    :param spacer: Kind of spacer to use ('n' for UnknownSeq spacer, 'stop' for all-frame stop codon spacer)

    :return: A single SeqRecord that is the product of the merge.
    """

    if spacer not in ('n', 'stop'):
        raise ValueError("Invalid spacer: %r, use either 'n' or 'stop'" % spacer)

    if not len(records):
        raise ValueError("No records given")

    if spacer == 'stop':
        spacer_seq = Seq(ALL_FRAME_STOP_MOTIF * 40 * length)
    else:
        spacer_seq = Seq('N' * 1000)

    new_rec = records[0]

    if len(records) == 1:
        return new_rec

    rec_id = new_rec.id
    rec_name = new_rec.name
    rec_desc = new_rec.description
    date = new_rec.annotations.get('date', '')
    source = new_rec.annotations.get("source", '')
    organism = new_rec.annotations.get('organism', '')
    taxonomy = new_rec.annotations.get('taxonomy', [])
    data_file_division = new_rec.annotations.get('data_file_division', 'UNK')
    topology = new_rec.annotations.get('topology', 'linear')

    for i, rec in enumerate(records[1:]):
        spacer_id = 'spacer_{}'.format(i + 1)

        spacer_feature = SeqFeature(FeatureLocation(0, int(length * 1000), 0),
                                    type='misc_feature', id=spacer_id,
                                    qualifiers={'note': [spacer_id]})

        spacer_rec = SeqRecord(spacer_seq, id=spacer_id, name=spacer_id,
                               description=spacer_id, features=[spacer_feature])

        new_rec = new_rec + spacer_rec + rec

    new_rec.id = rec_id
    new_rec.name = rec_name
    new_rec.description = rec_desc
    new_rec.annotations['molecule_type'] = "DNA"
    new_rec.annotations["date"] = date
    new_rec.annotations["source"] = source
    new_rec.annotations["organism"] = organism
    new_rec.annotations["taxonomy"] = taxonomy
    new_rec.annotations["data_file_division"] = data_file_division
    new_rec.annotations["topology"] = topology

    return new_rec



gene_list = sys.argv[1]
gene_names = ""
gb_file = ""

with open(gene_list,'r') as input_file:
    gene_names=[line.strip('\n') for line in input_file]

if len(sys.argv) > 2:
    for gb in sys.argv[2:]:
        print("arguments: ", gb)
        if gb.endswith(".gz"):
            #check if file is gzip compressed
            #handle = gzip.open(gb, "rt")
            with gzip.open(gb, 'rt') as gbin:
                with open('tmp_gb.gbff', 'wt') as gbout:
                    data = gbin.read()
                    gbout.write(data)
                    gb_file = "tmp_gb.gbff"
                    gbout.close()
        else:
            gb_file=gb
        #gb_object=SeqIO.read(gb_file,'gb')
        records = []
        #allgenes = []
        #records.extend(SeqIO.parse(seqfile, 'genbank'))
        #with open(gb_file) as handle:
            #record = SeqIO.read(handle, 'gb')
        records.extend(SeqIO.parse(gb_file, 'genbank'))
        merged_record = merge(records, 10, 'n')
        outfile = open('tmp_gb_merged.gbff', 'w')
        SeqIO.write([merged_record], outfile, 'genbank')
        outfile.close()
        gb_file = 'tmp_gb_merged.gbff'
        with open(gb_file) as handle:
            record = SeqIO.read(gb_file, 'genbank')
            allgenes=[feature for feature in record.features if feature.type =='gene']
            print("total number of genes:", len(allgenes))
        gene_sequences=[]
        # code inspired by the following GitHub page: https://github.com/vappiah/Python-Bioinformatics-Hacks/tree/main/Notebooks 
        gb_object=SeqIO.read(gb_file,'gb')    
        for gene in allgenes:
            if 'gene' in gene.qualifiers.keys():
                gene_name=gene.qualifiers['gene'][0]
                if gene_name in gene_names:
                    extract=gene.extract(gb_object)
                    extract.id=gene_name
                    extract.description=''
                    extract.annotation={"molecule_type":"DNA"}
                    gene_sequences.append(extract)
                    print('gene %s has been found'%gene_name)
        #len(gene_sequences)
        file_name = os.path.basename(gb)
        res_fasta = os.path.splitext(file_name)[0]+'_spec_genes.fasta'
        res_gb = os.path.splitext(file_name)[0]+'_spec_genes.gbff'
        result_fasta = open(res_fasta, 'wt')
        #result_gb = open(res_gb, 'wt')
        SeqIO.write(gene_sequences, result_fasta, 'fasta')
        #print("gene sequences = ", gene_sequences)
        #SeqIO.write(gene_sequences[0], result_gb, 'genbank')
        for g in gene_sequences:
            g.annotations={"molecule_type": "DNA"}
        recordseq = []
        with open (res_gb, "w") as outfile:
            #SeqIO.write([gene for gene in gene_sequences], outfile, "genbank")
            #recordseq.extend(SeqIO.parse(outfile, 'genbank'))
            #merged_recordseq = merge(recordseq, 20, 'n')
            SeqIO.write(gene_sequences, outfile, 'genbank')
            outfile.close()
        recordseq = []
        recordseq.extend(SeqIO.parse(res_gb, 'genbank'))
        merged_recordseq = merge(recordseq, 10, 'n')
        #outfile.close()
        outfile2 = open(res_gb, 'w')
        SeqIO.write([merged_recordseq], outfile2, 'genbank')
        outfile2.close()
        #os.remove('./tmp_gb.gbff')
        #os.remove('./tmp_gb_merged.gbff')

if os.path.exists('./tmp_gb.gbff'):
    os.remove('./tmp_gb.gbff')
if os.path.exists('./tmp_gb_merged.gbff'):
    os.remove('./tmp_gb_merged.gbff')
