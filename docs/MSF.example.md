# MAPI Sample Folders (MSFs): Minimal Examples

## Create / enter a sample folder

```bash
cd ~/project/
mkdir mapi-sampleA
cd mapi-sampleA
```

Put inputs in the folder (any filenames are fine):

```bash
ls
reads_R1.fastq.gz
reads_R2.fastq.gz
```

## Run a module (no paths, no output flags)

```bash
mapi modules bwa --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz
```

## What happens automatically

Because you are inside `mapi-sampleA/`, MAPI treats this as an MSF and uses a standard backend output folder:

```
mapi-sampleA/
├── reads_R1.fastq.gz
├── reads_R2.fastq.gz
└── sampleA-output/
    ├── sampleA.bwa_bam.bam
    ├── sampleA.bwa_bam.bam.bai
    └── sampleA.bwa_log.log
```

No output directory specified.  
No filenames invented by the user.  
No paths typed.

## Reusing outputs by *name*, not path (future-facing)

Later modules can refer to outputs logically (example):

```bash
mapi modules summarize_bam bwa_bam
```

Instead of:

```bash
mapi modules summarize_bam sampleA-output/sampleA.bwa_bam.bam
```
