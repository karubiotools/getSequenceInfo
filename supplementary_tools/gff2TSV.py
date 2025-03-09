import csv, os
import argparse

def parse_gff3(gff3_files, output_format='tsv'):
    for gff in gff3_files:
        basename = os.path.basename(gff)+'.tsv'
        with open(gff, 'r') as infile, open(basename, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t' if output_format == 'tsv' else ',')
            writer.writerow(["feature", "start", "stop", "strand", "gene_name"])
        
            for line in infile:
                if line.startswith('#'):
                    continue  # Skip comments and headers
            
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue  # Skip malformed lines
            
                feature, start, stop, strand, attributes = fields[2], fields[3], fields[4], fields[6], fields[8]
            
                gene_name = "HP"
                attributes_dict = {key: value for key, value in (attr.split('=') for attr in attributes.split(';') if '=' in attr)}
                if 'Name' in attributes_dict:
                    gene_name = attributes_dict['Name']
                elif 'gene' in attributes_dict:
                    gene_name = attributes_dict['gene']
                #elif 'ID' in attributes_dict:
                 #   gene_name = attributes_dict['ID']

                if feature == 'gene':
                    writer.writerow([feature, start, stop, strand, gene_name])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a GFF3 file to a TSV or CSV file.")
    parser.add_argument("-i", "--input_gff3", nargs='+', default="sequence.gff3", help="Input GFF3 file(s)", required=True)
    #parser.add_argument("--output_file", default="output.tsv", help="Output file (TSV or CSV)")
    parser.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="Output format (default: TSV)")
    args = parser.parse_args()
    
    parse_gff3(args.input_gff3, args.format)
