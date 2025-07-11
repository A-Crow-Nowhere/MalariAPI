import sys
import time
from collections import Counter, defaultdict
import json
from Bio import Entrez
import os

Entrez.email = "brownnoah456@gmail.com"  # << Replace with your email here!

def extract_accession(sseqid):
    """
    Extract accession from sseqid field like 'ref|NC_037282.1|'
    Returns 'NC_037282.1' or 'Unknown' if not parsable.
    """
    if sseqid in ["N/A", ".", "-", "", "Unknown"]:
        return "Unknown"
    if "|" in sseqid:
        parts = sseqid.split("|")
        if len(parts) > 1 and parts[1]:
            return parts[1]
        else:
            return parts[-1]  # fallback
    return sseqid.strip()

def fetch_genus_for_accessions(accessions):
    genuses = {}
    batch_size = 50
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        ids = ",".join(batch)
        print(f"  - Querying batch {i//batch_size + 1}: {len(batch)} accessions")
        try:
            handle = Entrez.esummary(db="nucleotide", id=ids, retmode="json")
            result = json.load(handle)
            handle.close()

            for acc in batch:
                if acc in result['result']:
                    organism = result['result'][acc].get('organism', '')
                    genus = organism.split()[0] if organism else "Unknown"
                    genuses[acc] = genus
                else:
                    print(f"    {acc} not found in Entrez result, assigned Unknown")
                    genuses[acc] = "Unknown"
        except Exception as e:
            print(f"⚠️ Failed to fetch genus for {ids}: {e}")
            for acc in batch:
                genuses[acc] = "Unknown"
    return genuses


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} blast_output_file", file=sys.stderr)
        sys.exit(1)

    blast_file = sys.argv[1]

    sample_name = os.path.basename(blast_file).replace('.nt.blast.out', '')

    print(f"⚙️ Parsing BLAST file '{blast_file}'...")
    acc_list = []
    read_to_acc = {}
    with open(blast_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            read_id = fields[0]
            sseqid = fields[1]
            acc = extract_accession(sseqid)
            if read_id not in read_to_acc:
                read_to_acc[read_id] = acc
                acc_list.append(acc)

    print(f"🔍 Found {len(read_to_acc)} unique reads.")

    unique_acc = list(set(acc_list))
    print(f"⚙️ Fetching genus info for {len(unique_acc)} unique accessions in batches of 50...")
    acc_to_genus = fetch_genus_for_accessions(unique_acc)

    genus_counts = Counter()
    unknown_genus_reads = 0
    for acc in acc_list:
        genus = acc_to_genus.get(acc, "Unknown")
        if genus == "Unknown":
            unknown_genus_reads += 1
        genus_counts[genus] += 1

    total_reads = len(acc_list)

    # Prepare TSV lines: header + data
    output_lines = ["genus\tpercentage"]
    for genus, count in genus_counts.most_common():
        percent = 100 * count / total_reads
        output_lines.append(f"{genus}\t{percent:.2f}")

    output_text = "\n".join(output_lines)
    print(output_text)

    # Write to file in blast_nt_out directory
    outdir = "blast_nt_out"
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, f"{sample_name}.genus_proportions.tsv")
    with open(output_file, "w") as outf:
        outf.write(output_text + "\n")

    print(f"\n✅ TSV output written to {output_file}")


if __name__ == "__main__":
    main()
