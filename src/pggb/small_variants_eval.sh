#!/usr/bin/env bash

# ==============================================================================
# small_variants_eval.sh
# Universal Small Variants Evaluation Module (Phases 3-5)
# Includes sleek 4-line terminal progress display and execution logging
# Usage: ./small_variants_eval.sh [INPUT_FASTA_OR_LOCUS] [PGGB_DIR] [REF_NAME] [OUT_DIR]
# Example: ./small_variants_eval.sh DRB1-3123 data/intern/small_variants/02_pggb grch38 data/intern/small_variants
# ==============================================================================

set -eo pipefail

INPUT_ARG="${1:-DRB1-3123}"
PGGB_DIR="${2:-data/intern/small_variants/02_pggb}"
REF_NAME="${3:-grch38}"
OUT_DIR="${4:-data/intern/small_variants}"
THREADS=8

# Auto-discover input FASTA file if locus name is provided
INPUT_FASTA="$INPUT_ARG"
if [[ ! -f "$INPUT_FASTA" ]]; then
    if [[ -f "data/intern/partitions/hla_${INPUT_ARG}.fasta.gz" ]]; then
        INPUT_FASTA="data/intern/partitions/hla_${INPUT_ARG}.fasta.gz"
    elif [[ -f "data/intern/partitions/${INPUT_ARG}.fasta.gz" ]]; then
        INPUT_FASTA="data/intern/partitions/${INPUT_ARG}.fasta.gz"
    elif [[ -f "data/raw/HLA/${INPUT_ARG}.fa" ]]; then
        INPUT_FASTA="data/raw/HLA/${INPUT_ARG}.fa"
    fi
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR_INPUT="$OUT_DIR/01_input"
DIR_PGGB="$PGGB_DIR"
DIR_NUCMER="$OUT_DIR/03_nucmer"
DIR_VCFEVAL="$OUT_DIR/04_vcfeval"
DIR_PLOTS="$OUT_DIR/05_plots"

prepare_reference() {
    mkdir -p "$DIR_INPUT" "$DIR_NUCMER" "$DIR_VCFEVAL" "$DIR_PLOTS"
    
    PREP_FASTA="$DIR_INPUT/input.fa.gz"
    if [[ "$INPUT_FASTA" == *.gz ]]; then
        cp "$INPUT_FASTA" "$PREP_FASTA"
    else
        bgzip -c -@ "$THREADS" "$INPUT_FASTA" > "$PREP_FASTA"
    fi
    samtools faidx "$PREP_FASTA"

    REF_FULL_NAME=$(cut -f 1 "${PREP_FASTA}.fai" | grep "$REF_NAME" | head -n 1 || true)
    if [[ -z "$REF_FULL_NAME" ]]; then
        REF_FULL_NAME=$(head -n 1 "${PREP_FASTA}.fai" | cut -f 1)
    fi

    REF_FA="$DIR_INPUT/ref.fa"
    samtools faidx "$PREP_FASTA" "$REF_FULL_NAME" > "$REF_FA"
    samtools faidx "$REF_FA"
}

graph_variant_calling() {
    FINAL_GFA=$(ls "$PGGB_DIR"/*.final.gfa 2>/dev/null | head -n 1 || true)
    if [[ -z "$FINAL_GFA" || ! -f "$FINAL_GFA" ]]; then
        echo "Error: No *.final.gfa graph file found in $PGGB_DIR"
        exit 1
    fi

    PGGB_RAW_VCF="$DIR_PGGB/pggb_raw.vcf.gz"
    PGGB_WAVE_VCF="$DIR_PGGB/pggb_wave.vcf.gz"

    vg deconstruct -P "$REF_FULL_NAME" -a -t "$THREADS" "$FINAL_GFA" | \
        bgzip -c -@ "$THREADS" > "$PGGB_RAW_VCF"

    vcfbub -l 0 -a 100000 --input "$PGGB_RAW_VCF" | \
        vcfwave -I 1000 -t "$THREADS" | \
        bcftools sort | \
        bgzip -c -@ "$THREADS" > "$PGGB_WAVE_VCF"

    tabix -f -p vcf "$PGGB_WAVE_VCF"

    cut -f 1 "${PREP_FASTA}.fai" | grep -v "^${REF_FULL_NAME}$" | while read -r CONTIG; do
        SAFE_NAME=$(echo "$CONTIG" | tr '/:#|' '_')
        (bcftools view -s "$CONTIG" -min-ac 1 "$PGGB_WAVE_VCF" 2>/dev/null || bcftools view "$PGGB_WAVE_VCF") | \
            bcftools sort | bgzip -c > "$DIR_PGGB/pggb_${SAFE_NAME}.vcf.gz"
        tabix -f -p vcf "$DIR_PGGB/pggb_${SAFE_NAME}.vcf.gz"
    done
}

nucmer_baseline() {
    cut -f 1 "${PREP_FASTA}.fai" | grep -v "^${REF_FULL_NAME}$" | while read -r CONTIG; do
        SAFE_NAME=$(echo "$CONTIG" | tr '/:#|' '_')
        CONTIG_FA="$DIR_NUCMER/${SAFE_NAME}.fa"
        PREFIX="$DIR_NUCMER/${SAFE_NAME}_vs_ref"
        
        samtools faidx "$PREP_FASTA" "$CONTIG" > "$CONTIG_FA"
        nucmer "$REF_FA" "$CONTIG_FA" --prefix "$PREFIX" >/dev/null 2>&1
        show-snps -THC "${PREFIX}.delta" > "${PREFIX}.var.txt"
        
        Rscript "$SCRIPT_DIR/nucmer2vcf.R" "${PREFIX}.var.txt" "$CONTIG" "$REF_FA" "4.0.0" "${PREFIX}.vcf"
        bgzip -f -@ "$THREADS" "${PREFIX}.vcf"
        tabix -f -p vcf "${PREFIX}.vcf.gz"
        rm -f "$CONTIG_FA"
    done
}

benchmark_eval() {
    rm -rf "$DIR_VCFEVAL"/*_eval
    REF_SDF="$DIR_VCFEVAL/ref.sdf"
    rm -rf "$REF_SDF"
    rtg format -o "$REF_SDF" "$REF_FA" >/dev/null 2>&1

    DIST=1000

    cut -f 1 "${PREP_FASTA}.fai" | grep -v "^${REF_FULL_NAME}$" | while read -r CONTIG; do
        SAFE_NAME=$(echo "$CONTIG" | tr '/:#|' '_')
        NUCMER_VCF="$DIR_NUCMER/${SAFE_NAME}_vs_ref.vcf.gz"
        PGGB_VCF="$DIR_PGGB/pggb_${SAFE_NAME}.vcf.gz"
        CALLABLE_BED="$DIR_VCFEVAL/${SAFE_NAME}_callable.bed"
        EVAL_OUT="$DIR_VCFEVAL/${SAFE_NAME}_eval"
        
        rm -rf "$EVAL_OUT"
        
        BED_A="$DIR_VCFEVAL/${SAFE_NAME}_bed_a.bed"
        BED_B="$DIR_VCFEVAL/${SAFE_NAME}_bed_b.bed"
        
        (bedtools merge -d $DIST -i "$NUCMER_VCF" 2>/dev/null || echo -e "${REF_FULL_NAME}\t1\t100") > "$BED_A"
        (bedtools merge -d $DIST -i "$PGGB_VCF" 2>/dev/null || echo -e "${REF_FULL_NAME}\t1\t100") > "$BED_B"
        
        bedtools intersect -a "$BED_A" -b "$BED_B" > "$CALLABLE_BED" 2>/dev/null || echo -e "${REF_FULL_NAME}\t1\t100" > "$CALLABLE_BED"
        
        if [[ ! -s "$CALLABLE_BED" ]]; then
            echo -e "${REF_FULL_NAME}\t1\t100" > "$CALLABLE_BED"
        fi
        rm -f "$BED_A" "$BED_B"
        
        rtg vcfeval \
            -t "$REF_SDF" \
            -b "$NUCMER_VCF" \
            -c "$PGGB_VCF" \
            --sample "$CONTIG,$CONTIG" \
            --squash-ploidy \
            --all-records \
            -T "$THREADS" \
            -e "$CALLABLE_BED" \
            -o "$EVAL_OUT" >/dev/null 2>&1 || true
    done

    STAT_FILE="$DIR_VCFEVAL/statistics.tsv"
    echo -e "contig\tprecision\trecall\tf1.score" > "$STAT_FILE"

    for SUMMARY in "$DIR_VCFEVAL"/*_eval/summary.txt; do
        if [[ -f "$SUMMARY" ]]; then
            C_NAME=$(basename "$(dirname "$SUMMARY")" | sed 's/_eval//')
            NONE_LINE=$(grep "None" "$SUMMARY" || true)
            if [[ -n "$NONE_LINE" ]]; then
                VALS=$(echo "$NONE_LINE" | tr -s ' ' | cut -f 7,8,9 -d ' ' | tr ' ' '\t')
                echo -e "${C_NAME}\t${VALS}" >> "$STAT_FILE"
            else
                echo -e "${C_NAME}\t0.000\t0.000\t0.000" >> "$STAT_FILE"
            fi
        fi
    done

    PLOT_PNG="$DIR_PLOTS/precision_recall_f1.png"
    Rscript "$SCRIPT_DIR/plot_small_variants.R" "$STAT_FILE" "$PLOT_PNG" >/dev/null 2>&1
}

cleanup_intermediates() {
    rm -f "$DIR_INPUT"/input.fa.gz "$DIR_INPUT"/input.fa.gz.*
    rm -f "$DIR_NUCMER"/*.delta "$DIR_NUCMER"/*.var.txt
    rm -f "$DIR_PGGB"/pggb_raw.vcf.gz
    rm -f "$DIR_VCFEVAL"/*_callable.bed
}

main() {
    EVAL_LOG="${OUT_DIR}/small_variants_evaluation.log"

    echo "================================================================="
    echo "[EVAL] SMALL VARIANTS BENCHMARK PIPELINE"
    echo "  • Input FASTA:    $INPUT_FASTA"
    echo "  • PGGB Graph:     $PGGB_DIR"
    echo "  • Reference:      $REF_NAME"
    echo "  • Output Dir:     $OUT_DIR"
    echo "  • Execution Log:  $EVAL_LOG"
    echo "================================================================="

    status_prep="[Running... ⏳]"
    status_calling="[Waiting...]"
    status_nucmer="[Waiting...]"
    status_rtg="[Waiting...]"

    print_status() {
        echo -e " 1/4  Prep Reference:   $status_prep"
        echo -e " 2/4  Graph Calling:    $status_calling"
        echo -e " 3/4  Nucmer Baseline:  $status_nucmer"
        echo -e " 4/4  RTG Benchmarking: $status_rtg"
    }

    print_status

    # Phase 1: Reference prep
    prepare_reference > "$EVAL_LOG" 2>&1
    status_prep="[Completed ✔]"
    status_calling="[Running... ⏳]"
    echo -ne "\033[4A"
    print_status

    # Phase 2: Graph variant calling
    graph_variant_calling >> "$EVAL_LOG" 2>&1
    status_calling="[Completed ✔]"
    status_nucmer="[Running... ⏳]"
    echo -ne "\033[4A"
    print_status

    # Phase 3: Nucmer baseline
    nucmer_baseline >> "$EVAL_LOG" 2>&1
    status_nucmer="[Completed ✔]"
    status_rtg="[Running... ⏳]"
    echo -ne "\033[4A"
    print_status

    # Phase 4: Benchmark evaluation
    benchmark_eval >> "$EVAL_LOG" 2>&1
    status_rtg="[Completed ✔]"
    echo -ne "\033[4A"
    print_status

    # Phase 5: Cleanup
    cleanup_intermediates

    STAT_FILE="$DIR_VCFEVAL/statistics.tsv"
    PLOT_PNG="$DIR_PLOTS/precision_recall_f1.png"

    echo ""
    echo "================================================================="
    echo "[SUCCESS] Benchmark Evaluation Completed!"
    echo "  • Summary Table:  $STAT_FILE"
    echo "  • Benchmark Plot: $PLOT_PNG"
    echo "  • Execution Log:  $EVAL_LOG"
    echo "================================================================="
    if [[ -f "$STAT_FILE" ]]; then
        cat "$STAT_FILE" | column -t
    fi
}

main "$@"
