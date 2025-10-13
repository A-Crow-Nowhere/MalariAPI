#!/usr/bin/env python3
import sys
import os
import time
import json
from collections import Counter
from Bio import Entrez

Entrez.email = "EMAIL"  # Replace this with your email

def extract_accession(sseqid):
    if sseqid in ["N/A", ".", "-", "", "Unknown"]:
        return "Unknown"
    if "|" in sseqid:
        parts = sseqid.split("|")
        return parts[1] if len(parts) > 1 and parts[1] else parts[-1]
    return sseqid.strip()

def fetch_genus_for_accessions(accessions):
    genuses = {}
    batch_size = 50
    print(f"⚙️ Fetching genus info for {len(accessions)} unique accessions in batches of {batch_size}...")

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        accession_to_uid = {}

        # Step 1: Accessions → UIDs
        for acc in batch:
            try:
                search_handle = Entrez.esearch(db="nucleotide", term=f"{acc}[Accession]", retmode="json")
                search_results = json.load(search_handle)
                search_handle.close()
                uid_list = search_results.get("esearchresult", {}).get("idlist", [])
                if uid_list:
                    accession_to_uid[acc] = uid_list[0]
                else:
                    print(f"    ❌ Accession {acc} could not be resolved to UID")
                    genuses[acc] = "Unknown"
            except Exception as e:
                print(f"    ⚠️ Failed to search UID for {acc}: {e}")
                genuses[acc] = "Unknown"

        uids = list(accession_to_uid.values())
        if not uids:
            continue

        # Step 2: UID → Genus
        try:
            summary_handle = Entrez.esummary(db="nucleotide", id=",".join(uids), retmode="json")
            summary_results = json.load(summary_handle)
            summary_handle.close()
            result_data = summary_results.get("result", {})

            for acc, uid in accession_to_uid.items():
                entry = result_data.get(uid, {})
                organism = entry.get("organism", "")
                genus = organism.split()[0] if organism else "Unknown"
                genuses[acc] = genus
                if genus == "Unknown":
                    print(f"    ⚠️ Genus not found for {acc}")
        except Exception as e:
            print(f"⚠️ Failed to fetch genus for UIDs: {uids}: {e}")
            for acc in accession_to_uid.keys():
                genuses[acc] = "Unknown"

        time.sleep(2.04)  # NCBI rate limit

    return genuses

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <blast_output_file>", file=sys.stderr)
        sys.exit(1)

    blast_file = sys.argv[1]
    sample_name = os.path.basename(blast_file).replace(".nt.blast.out", "")
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
            read_id, sseqid = fields[0], fields[1]
            acc = extract_accession(sseqid)
            if read_id not in read_to_acc:
                read_to_acc[read_id] = acc
                acc_list.append(acc)

    print(f"🔍 Found {len(read_to_acc)} unique reads.")

    total_reads = len(acc_list)
    acc_to_genus = fetch_genus_for_accessions(sorted(set(acc_list)))

    genus_counts = Counter()
    for acc in acc_list:
        genus = acc_to_genus.get(acc, "Unknown")
        genus_counts[genus] += 1

    # Create output table
    output_lines = ["genus\tpercentage\tnumber"]
    for genus, count in genus_counts.most_common():
        percent = 100 * count / total_reads
        output_lines.append(f"{genus}\t{percent:.2f}\t{count}")

    # Save output
    os.makedirs("blastOut", exist_ok=True)
    output_file = os.path.join("blastOut", f"{sample_name}.genus_proportions.tsv")
    with open(output_file, "w") as outf:
        outf.write("\n".join(output_lines) + "\n")

    print(f"✅ TSV output written to {output_file}")

if __name__ == "__main__":
    main()
