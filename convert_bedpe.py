"""Convert 10x Genomics feature_linkage.bedpe to peak_gene_links.csv format"""

import re

# Input/Output paths
BEDPE_PATH = "2_data_intermediate/atac_reference/analysis/feature_linkage/feature_linkage.bedpe"
CSV_PATH = "2_data_intermediate/atac_reference/peak_gene_links.csv"

# Parse BEDPE and convert
print("Parsing BEDPE file...")
links = []

with open(BEDPE_PATH, 'r') as f:
    for i, line in enumerate(f):
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        
        # BEDPE format: chr1, start1, end1, chr2, start2, end2, name, score, ...
        chr1, start1, end1 = parts[0], parts[1], parts[2]
        chr2, start2, end2 = parts[3], parts[4], parts[5]
        name = parts[6]  # Format: <GENE_promoter><GENE> or similar
        score = float(parts[7]) if parts[7] else 0.0
        
        # Extract gene name from the name field
        gene_match = re.search(r'<([^>]+)>', name)
        gene = gene_match.group(1).replace('_promoter', '') if gene_match else 'Unknown'
        
        # Create peak identifier
        peak = f"{chr1}:{start1}-{end1}"
        
        # Calculate distance
        distance = abs(int(start2) - int(start1))
        
        # Cell type is "PBMC" (reference data)
        cell_type = "PBMC"
        
        links.append({
            'peak': peak,
            'gene': gene,
            'cell_type': cell_type,
            'distance_bp': distance,
            'score': score
        })
        
        if i > 0 and i % 100000 == 0:
            print(f"  Processed {i:,} lines...")

print(f"Total links parsed: {len(links):,}")

# Write CSV
print(f"Writing to {CSV_PATH}...")
with open(CSV_PATH, 'w') as f:
    f.write("peak,gene,cell_type,distance_bp,score\n")
    for link in links:
        f.write(f"{link['peak']},{link['gene']},{link['cell_type']},{link['distance_bp']},{link['score']}\n")

print(f"Done! Created {CSV_PATH} with {len(links):,} peak-gene links")
