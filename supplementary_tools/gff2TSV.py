import csv
import argparse

def parse_gff3(gff3_file, output_file, output_format='tsv'):
    with open(gff3_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t' if output_format == 'tsv' else ',')
        writer.writerow(["feature", "start", "stop", "strand", "gene_name"])
        
        for line in infile:
            if line.startswith('#'):
                continue  # Skip comments and headers
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # Skip malformed lines
            
            feature, start, stop, strand, attributes = fields[2], fields[3], fields[4], fields[6], fields[8]
            
            gene_name = "NA"
            attributes_dict = {key: value for key, value in (attr.split('=') for attr in attributes.split(';') if '=' in attr)}
            if 'Name' in attributes_dict:
                gene_name = attributes_dict['Name']
            elif 'gene' in attributes_dict:
                gene_name = attributes_dict['gene']
            elif 'ID' in attributes_dict:
                gene_name = attributes_dict['ID']
            
            writer.writerow([feature, start, stop, strand, gene_name])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a GFF3 file to a TSV or CSV file.")
    parser.add_argument("--input_gff3", default="sequence.gff3", help="Input GFF3 file")
    parser.add_argument("--output_file", default="output.tsv", help="Output file (TSV or CSV)")
    parser.add_argument("--format", choices=["tsv", "csv"], default="tsv", help="Output format (default: TSV)")
    args = parser.parse_args()
    
    parse_gff3(args.input_gff3, args.output_file, args.format)
