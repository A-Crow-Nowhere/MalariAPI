## FOR MODULE WRAPPING:

### Strictness & Safety
•	set -euo pipefail is mandatory.
•	Quote all paths.
•	Do not cd; use absolute/derived paths.
•	Fail fast with clear messages; let the pipeline capture logs.
### Standard Inputs & Outputs
The pipeline (mapi run) assumes only these conventions:
#### Reads Filtering Modules
•	Inputs: -1 <R1.fq.gz> -2 <R2.fq.gz> or -s <SE.fq.gz>
•	Outputs:
o	Paired:
${OUTDIR}/${PREFIX}_R1.filtered.fastq.gz
${OUTDIR}/${PREFIX}_R2.filtered.fastq.gz
o	Single:
${OUTDIR}/${PREFIX}.filtered.fastq.gz
•	Reports (if any): ${OUTDIR}/${PREFIX}.<tool>.{json,html}
#### QC Modules (e.g., FastQC)
•	Inputs: 1+ FASTQ(.gz)
•	Outputs: native tool output into ${OUTDIR}/ (no rename required)
#### Alignment Modules (e.g., BWA-MEM2)
•	Inputs: filtered reads + -r|--ref <ref.fa>
•	Outputs:
o	${OUTDIR}/${PREFIX}.sorted.bam
o	${OUTDIR}/${PREFIX}.sorted.bam.bai
•	Behavior:
o	Auto-index reference if needed.
o	Pipe to samtools sort; then samtools index.
### Variant / SV Modules
•	SNP/indel: ${OUTDIR}/${PREFIX}.vcf.gz (+ .tbi)
•	SV: ${OUTDIR}/${PREFIX}.sv.vcf.gz (+ .tbi)
•	Sidecars allowed, keep ${PREFIX}.* pattern.
### Logging & Idempotency
•	Keep stdout minimal (status lines). The pipeline redirects to ./<prefix>_MAPI/logs/NN_<step>.log.
•	Do not overwrite existing primary artifacts unless you expose --force.
•	Create ${OUTDIR}/_mapi_done on success (optional sentinel).
### Resource Controls
•	Honor THREADS (from env); provide a sane default when not set.
•	If you allocate temp files, use ${OUTDIR} or mktemp -d.
Reserved Artifact Labels (for documentation & future manifesting)
•	reads/raw, reads/filtered, qc/report
•	alignment/bam
•	variants/vcf, sv/vcf
•	metrics/json

