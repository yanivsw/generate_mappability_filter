# This script is used to generate a mappability filter for a given fasta file using Heng Li's
# seqbility tools (https://lh3lh3.users.sourceforge.net/snpable.shtml) to generate a mappability filter.
# The script uses the following arguments:
# -i, --input_fasta
# -k, --kmer_length
# -r, --stringency (i.e. the number of mismatches allowed)
# The script requires bwa, samtools and parallel to be installed on the system.

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 -i=<input_fasta> -k=<kmer_length> -r=<stringency>"
    echo
    echo "Arguments:"
    echo "  -i, --input_fasta    Path to the input FASTA file."
    echo "  -k, --kmer_length    Length of the k-mers to use."
    echo "  -r, --stringency     Number of mismatches allowed."
    echo
    echo "Description:"
    echo "  This script generates a mappability filter for a given FASTA file."
    echo "  It requires bwa, samtools, and parallel to be installed on the system."
    exit 0
fi

for arg in "$@"; do
    case $arg in
        -i=*|--input_fasta=*)
            input_fasta="${arg#*=}"
            ;;
        -k=*|--kmer_length=*)
            kmer_length="${arg#*=}"
            ;;
        -r=*|--stringency=*)
            stringency="${arg#*=}"
            ;;
        *)
            echo "Unknown option $arg"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$input_fasta" || -z "$kmer_length" || -z "$stringency" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 -i=<input_fasta> -k=<kmer_length> -r=<stringency>"
    exit 1
fi

echo $(dirname $input_fasta)
echo $(basename $input_fasta)

output_fasta=${input_fasta/.fa/".mask_${kmer_length}_${stringency#*.}.fa"}
output_bed=${input_fasta/.fa/".mask_${kmer_length}_${stringency#*.}.bed"}
output_bam=${input_fasta/.fa/"_${kmer_length}.bam"}
tmp_dir=$(dirname $input_fasta)/tmp_mask_${kmer_length}_${stringency#*.}

echo $output_fasta
echo $output_bed
echo $output_bam
echo $tmp_dir

echo $kmer_length
echo $stringency

mkdir -p $tmp_dir

# Ensure BWA index exists
if [[ ! -f "${input_fasta}.bwt" || ! -f "${input_fasta}.pac" || ! -f "${input_fasta}.ann" || ! -f "${input_fasta}.amb" || ! -f "${input_fasta}.sa" ]]; then
    echo "BWA index files not found. Generating BWA index..."
    bwa index $input_fasta
fi

$(dirname "$0")/seqbility-20091110/splitfa $input_fasta $kmer_length | split -l 20000000 - $tmp_dir/x

parallel "bwa aln -R 1000000 -O 3 -E 3 $input_fasta {} > {.}.sai" ::: $tmp_dir/x*

parallel "bwa samse $input_fasta {} {.} > {.}.sam" ::: $tmp_dir/x*.sai

parallel "samtools view -bS {} > {.}.bam" ::: $tmp_dir/*.sam

samtools cat -o $output_bam $tmp_dir/x*.bam

samtools view $output_bam | $(dirname "$0")/seqbility-20091110/gen_raw_mask.pl > $(dirname $input_fasta)/rawMask_${kmer_length}.fa

$(dirname "$0")/seqbility-20091110/gen_mask -l $kmer_length -r $stringency $(dirname $input_fasta)/rawMask_${kmer_length}.fa > $(dirname $input_fasta)/mask_${kmer_length}_${stringency#*.}.fa

$(dirname "$0")/seqbility-20091110/apply_mask_s $(dirname $input_fasta)/mask_${kmer_length}_${stringency#*.}.fa $input_fasta > $output_fasta

python3 $(dirname "$0")/convert_seqbility_fa_to_bed.py $output_fasta $output_bed
