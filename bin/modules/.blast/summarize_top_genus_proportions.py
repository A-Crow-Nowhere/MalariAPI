#!/usr/bin/env python3
import sys, os, time, json, argparse
from collections import Counter

def extract_accession(sseqid: str) -> str:
    if sseqid in ["N/A", ".", "-", "", "Unknown"]:
        return "Unknown"
    if "|" in sseqid:
        parts = sseqid.split("|")
        return parts[1] if len(parts) > 1 and parts[1] else parts[-1]
    return sseqid.strip()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blast_output_file")
    ap.add_argument("--entrez-email", default="")
    ap.add_argument("--no-entrez", action="store_true",
                    help="Do NOT query NCBI. Use sscinames column from BLAST output to infer genus.")
    args = ap.parse_args()

    blast_file = args.blast_output_file
    sample_name = os.path.basename(blast_file).replace(".nt.blast.out", "")
    print(f"⚙️ Parsing BLAST file '{blast_file}'...")

    # Read BLAST output: outfmt '6 qseqid sseqid pident length mismatch evalue sscinames'
    # columns: 0=qseqid 1=sseqid ... 5=evalue 6=sscinames
    read_to_acc = {}
    acc_list = []
    ssciname_list = []

    with open(blast_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            read_id = fields[0]
            if read_id in read_to_acc:
                continue

            sseqid = fields[1]
            acc = extract_accession(sseqid)
            read_to_acc[read_id] = acc
            acc_list.append(acc)

            ssciname = fields[6] if len(fields) > 6 else ""
            ssciname_list.append(ssciname)

    print(f"🔍 Found {len(read_to_acc)} unique reads.")

    total_reads = len(acc_list)
    genus_counts = Counter()

    if args.no_entrez:
        for title in ssciname_list:
            title = (title or "").strip()
            if not title or title.upper() == "N/A":
                genus = "Unknown"
            else:
                genus = title.split()[0]
            genus_counts[genus] += 1


    else:
        # Original behavior: Entrez lookups by accession
        try:
            from Bio import Entrez
        except Exception as e:
            print(f"❌ Biopython not available: {e}", file=sys.stderr)
            sys.exit(2)

        if args.entrez_email:
            Entrez.email = args.entrez_email
        else:
            # keep your original placeholder behavior, but warn
            Entrez.email = "EMAIL"
            print("⚠️ Entrez email not set. Pass --entrez-email you@domain for best practice.", file=sys.stderr)

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
                        h = Entrez.esearch(db="nucleotide", term=f"{acc}[Accession]", retmode="json")
                        search_results = json.load(h)
                        h.close()
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
                    h = Entrez.esummary(db="nucleotide", id=",".join(uids), retmode="json")
                    summary_results = json.load(h)
                    h.close()
                    result_data = summary_results.get("result", {})

                    for acc, uid in accession_to_uid.items():
                        entry = result_data.get(uid, {})
                        organism = entry.get("organism", "")
                        genus = organism.split()[0] if organism else "Unknown"
                        genuses[acc] = genus
                except Exception as e:
                    print(f"⚠️ Failed to fetch genus for UIDs: {uids}: {e}")
                    for acc in accession_to_uid.keys():
                        genuses[acc] = "Unknown"

                time.sleep(2.04)  # NCBI rate limit
            return genuses

        acc_to_genus = fetch_genus_for_accessions(sorted(set(acc_list)))
        for acc in acc_list:
            genus_counts[acc_to_genus.get(acc, "Unknown")] += 1

    # Output TSV (same columns/format as your original)
    output_lines = ["genus\tpercentage\tnumber"]
    for genus, count in genus_counts.most_common():
        percent = 100 * count / total_reads if total_reads else 0.0
        output_lines.append(f"{genus}\t{percent:.2f}\t{count}")

    # Write into same folder as blast output (MAPI-friendly, no implicit blastOut)
    out_dir = os.path.dirname(os.path.abspath(blast_file))
    output_file = os.path.join(out_dir, f"{sample_name}.genus_proportions.tsv")
    with open(output_file, "w") as outf:
        outf.write("\n".join(output_lines) + "\n")

    print(f"✅ TSV output written to {output_file}")

if __name__ == "__main__":
    main()
