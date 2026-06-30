import pyensembl
from collections import defaultdict

def main():
    ensembl = pyensembl.EnsemblRelease(release=111, species="human")
    ensembl.download()
    ensembl.index()
    
    conn = ensembl.db.connection
    cursor = conn.cursor()
    cursor.execute("SELECT gene_id, gene_name, start, end FROM exon;")
    exons = cursor.fetchall()
    
    gene_exons = defaultdict(list)
    gene_names = {}
    for gene_id, gene_name, start, end in exons:
        gene_exons[gene_id].append((start, end))
        gene_names[gene_id] = gene_name
        
    gene_lengths = {}
    for gene_id, intervals in gene_exons.items():
        intervals.sort(key=lambda x: x[0])
        merged = []
        for start, end in intervals:
            if not merged or merged[-1][1] < start:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        gene_lengths[gene_id] = sum(end - start + 1 for start, end in merged)
        
    # Check specific genes
    target_names = ["B2M", "CD28", "MIF", "HLA-A", "HLA-B", "HLA-C", "COL1A1"]
    name_to_id = {v: k for k, v in gene_names.items() if v in target_names}
    
    print("Genomic lengths of target genes:")
    for name in target_names:
        gid = name_to_id.get(name)
        if gid:
            print(f"Gene name: {name}, ID: {gid}, Union Exon Length: {gene_lengths[gid]} bp")
        else:
            # Let's search by name in all keys
            matches = [k for k, v in gene_names.items() if v == name]
            if matches:
                for gid in matches:
                    print(f"Gene name: {name}, ID: {gid}, Union Exon Length: {gene_lengths[gid]} bp")
            else:
                print(f"Gene name: {name} not found in exon DB.")

if __name__ == "__main__":
    main()
