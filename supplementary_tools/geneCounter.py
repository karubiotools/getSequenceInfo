import argparse
import gzip
import os
from collections import defaultdict
from Bio import SeqIO
import csv

def parse_genbank_file(file_path):
    """Parse un fichier GenBank et retourne un dictionnaire des gènes et leur occurrence."""
    genes = defaultdict(int)
    try:
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as handle:
                for record in SeqIO.parse(handle, 'genbank'):
                    for feature in record.features:
                        if feature.type == 'gene':
                            gene_name = feature.qualifiers.get('gene', ['Unknown'])[0]
                            genes[gene_name] += 1
        else:
            with open(file_path, 'r') as handle:
                for record in SeqIO.parse(handle, 'genbank'):
                    for feature in record.features:
                        if feature.type == 'gene':
                            gene_name = feature.qualifiers.get('gene', ['Unknown'])[0]
                            genes[gene_name] += 1
    except Exception as e:
        print(f"Erreur lors de la lecture du fichier {file_path}: {e}")
    return genes

def main():
    parser = argparse.ArgumentParser(description="Calculer l'occurrence des gènes dans des fichiers GenBank
.")
    parser.add_argument('files', nargs='+', help="Chemins des fichiers GenBank (.gbk, .gbff, .gb, .gz)")
    parser.add_argument('--genes', nargs='*', default=None, help="Liste des gènes à rechercher (optionnel)"
)
    parser.add_argument('--output', default='gene_occurrence.csv', help="Nom du fichier CSV de sortie (par
défaut: gene_occurrence.csv)")
    args = parser.parse_args()

    # Dictionnaire pour stocker les résultats
    results = defaultdict(dict)

    # Liste de tous les gènes uniques si aucune liste n'est fournie
    all_genes = set()

    # Lire chaque fichier GenBank
    for file_path in args.files:
        file_name = os.path.basename(file_path)
        genes = parse_genbank_file(file_path)
        results[file_name] = genes
        if args.genes is None:
            all_genes.update(genes.keys())

    # Si une liste de gènes est fournie, utiliser celle-ci
    if args.genes is not None:
        all_genes = set(args.genes)

    # Générer le tableau des occurrences
    table = []
    header = ['Genome'] + list(all_genes)
    table.append(header)

    for file_name, genes in results.items():
        row = [file_name]
        for gene in all_genes:
            row.append(genes.get(gene, 0))
        table.append(row)

    # Afficher le tableau
    for row in table:
        print("\t".join(map(str, row)))

    # Enregistrer le tableau dans un fichier CSV
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(table)

    print(f"Les résultats ont été enregistrés dans {args.output}")

if __name__ == "__main__":
    main()
