#!/usr/bin/env bash
set -euo pipefail

out_dir="${STAR_RS_YEAST_DIR:-.tmp/yeast_conformance}"
mkdir -p "${out_dir}"

fasta_gz="${out_dir}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
gtf_gz="${out_dir}/Saccharomyces_cerevisiae.R64-1-1.62.gtf.gz"
fasta="${fasta_gz%.gz}"
gtf="${gtf_gz%.gz}"
fastq="${out_dir}/SRR10143877.fastq"
spliced_fastq="${out_dir}/SRR10143877.spliced.fastq"
subset_chrom="${STAR_RS_YEAST_SUBSET_CHROM:-I}"
subset_fasta="${out_dir}/Saccharomyces_cerevisiae.R64-1-1.chrom_${subset_chrom}.fa"
subset_gtf="${out_dir}/Saccharomyces_cerevisiae.R64-1-1.chrom_${subset_chrom}.gtf"
star_bin="${STAR_RS_ORIGINAL_STAR:-STAR/bin/Linux_x86_64_static/STAR}"

fasta_url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-62/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
gtf_url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-62/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.62.gtf.gz"

fetch() {
  local url="$1"
  local output="$2"
  curl --ipv4 --fail --location --retry 4 --retry-delay 2 --connect-timeout 20 --max-time 300 \
    "${url}" -o "${output}"
}

if [[ ! -s "${fasta_gz}" ]]; then
  fetch "${fasta_url}" "${fasta_gz}"
fi
if [[ ! -s "${gtf_gz}" ]]; then
  fetch "${gtf_url}" "${gtf_gz}"
fi

if [[ ! -s "${fasta}" ]]; then
  gzip -dc "${fasta_gz}" > "${fasta}"
fi
if [[ ! -s "${gtf}" ]]; then
  gzip -dc "${gtf_gz}" > "${gtf}"
fi

if [[ ! -s "${fastq}" ]]; then
  fasterq-dump SRR10143877 --outdir "${out_dir}" --threads "${STAR_RS_YEAST_FASTQ_THREADS:-2}"
fi

if [[ ! -s "${subset_fasta}" || ! -s "${subset_gtf}" ]]; then
  python3 tools/make_yeast_subset.py "${fasta}" "${gtf}" "${subset_fasta}" "${subset_gtf}" "${subset_chrom}"
fi

if [[ "${STAR_RS_YEAST_SPLICE_ENRICHED:-1}" == "1" && ! -s "${spliced_fastq}" ]]; then
  if [[ ! -x "${star_bin}" ]]; then
    printf 'Skipping splice-enriched FASTQ creation because %s is not executable\n' "${star_bin}" >&2
  else
    splice_genome_dir="${out_dir}/splice_select_genome"
    splice_prefix="${out_dir}/splice_select/"
    mkdir -p "${splice_genome_dir}" "${splice_prefix}"
    if [[ ! -s "${splice_genome_dir}/Genome" ]]; then
      "${star_bin}" \
        --runMode genomeGenerate \
        --genomeDir "${splice_genome_dir}" \
        --genomeFastaFiles "${fasta}" \
        --sjdbGTFfile "${gtf}" \
        --sjdbOverhang "${STAR_RS_REAL_SJDB_OVERHANG:-50}" \
        --genomeSAindexNbases "${STAR_RS_SPLICE_SELECT_GENOME_SA_INDEX_NBASES:-5}" \
        --genomeChrBinNbits "${STAR_RS_REAL_GENOME_CHR_BIN_NBITS:-8}" \
        --limitGenomeGenerateRAM "${STAR_RS_REAL_LIMIT_GENOME_RAM:-1000000000}"
    fi
    "${star_bin}" \
      --genomeDir "${splice_genome_dir}" \
      --readFilesIn "${fastq}" \
      --readMapNumber "${STAR_RS_SPLICE_SELECT_SCAN_READS:-200000}" \
      --outSAMtype SAM \
      --outFileNamePrefix "${splice_prefix}"
    python3 tools/extract_fastq_by_sam_cigar.py \
      "${splice_prefix}/Aligned.out.sam" \
      "${fastq}" \
      "${spliced_fastq}" \
      "${STAR_RS_SPLICE_SELECT_READS:-100}"
  fi
fi

if [[ "${STAR_RS_YEAST_SPLICE_ENRICHED:-1}" == "1" && -s "${spliced_fastq}" ]]; then
  env_fasta="${fasta}"
  env_gtf="${gtf}"
  env_fastq="${spliced_fastq}"
  env_sa_nbases="\${STAR_RS_REAL_GENOME_SA_INDEX_NBASES:-5}"
  env_read_limit="\${STAR_RS_REAL_READ_LIMIT:-100}"
  env_min_splices="\${STAR_RS_REAL_MIN_SPLICES:-1}"
elif [[ "${STAR_RS_YEAST_FULL:-0}" == "1" ]]; then
  env_fasta="${fasta}"
  env_gtf="${gtf}"
  env_fastq="${fastq}"
  env_sa_nbases="\${STAR_RS_REAL_GENOME_SA_INDEX_NBASES:-5}"
  env_read_limit="\${STAR_RS_REAL_READ_LIMIT:-1000}"
  env_min_splices="\${STAR_RS_REAL_MIN_SPLICES:-1}"
else
  env_fasta="${subset_fasta}"
  env_gtf="${subset_gtf}"
  env_fastq="${fastq}"
  env_sa_nbases="\${STAR_RS_REAL_GENOME_SA_INDEX_NBASES:-3}"
  env_read_limit="\${STAR_RS_REAL_READ_LIMIT:-1000}"
  env_min_splices="\${STAR_RS_REAL_MIN_SPLICES:-0}"
fi

cat > "${out_dir}/env.sh" <<EOF
export STAR_RS_REAL_FASTA="${PWD}/${env_fasta}"
export STAR_RS_REAL_GTF="${PWD}/${env_gtf}"
export STAR_RS_REAL_FASTQ="${PWD}/${env_fastq}"
export STAR_RS_REAL_READ_LIMIT="${env_read_limit}"
export STAR_RS_REAL_GENOME_SA_INDEX_NBASES="${env_sa_nbases}"
export STAR_RS_REAL_GENOME_CHR_BIN_NBITS="\${STAR_RS_REAL_GENOME_CHR_BIN_NBITS:-8}"
export STAR_RS_REAL_LIMIT_GENOME_RAM="\${STAR_RS_REAL_LIMIT_GENOME_RAM:-1000000000}"
export STAR_RS_REAL_SJDB_OVERHANG="\${STAR_RS_REAL_SJDB_OVERHANG:-50}"
export STAR_RS_REAL_MIN_SPLICES="${env_min_splices}"
EOF

printf 'Prepared yeast conformance data in %s\n' "${out_dir}"
printf 'Using FASTA: %s\n' "${env_fasta}"
printf 'Using GTF:   %s\n' "${env_gtf}"
printf 'Using FASTQ: %s\n' "${env_fastq}"
printf 'Run:\n'
printf '  source %s/env.sh\n' "${out_dir}"
printf '  cargo test real_world_conformance_from_env_matches_original_star_core_sam_fields --test cpp_parity -- --nocapture\n'
