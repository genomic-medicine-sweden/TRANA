# Quit script if any of the commands fail. Note that && commands should be in parentheses for this to work.
set -eo pipefail

# Trap any errors and exit with a useful error message
trap 'exit_status="$?" && echo "Error occurred on line $LINENO: $BASH_COMMAND" && echo "Exit status: $exit_status" && exit "$exit_status"' ERR

# Usage message
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi


# Set the directory to the first command line argument and extract absolute path
directory="$1"
directory=$(realpath "$directory")

# Make an input samplesheet
#mkdir -p sample_sheets
ls -1 "$directory"/*fastq.gz > "reads.txt"
echo "sample,fastq_1,fastq_2" > "samplesheet_merged.csv"
while IFS= read -r line; do
    # Get the filename
    read1_n=$(basename "$line")
    # Full path to forward read
    read1_n_path="${directory}/${read1_n}"
    # Entry name (everything in the filename before ".fastq.gz")
    sample_n=$(echo "$read1_n" | sed  's/\.fastq\.gz//')
    # Append to the sample_sheet.csv file
    echo  "${sample_n},${read1_n_path}," >> "samplesheet_merged.csv"
    echo
done < "reads.txt"

cat "./samplesheet_merged.csv"
echo
echo "file saved as samplesheet_merged.csv"
echo
rm -f ./reads.txt
echo "Script finished"

