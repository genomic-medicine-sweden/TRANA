# quit script if any of the commands fails. Note that && commands should be in paranthesis for this to work.
set -eo pipefail
trap 'exit_status="$?" && echo Failed on line: $LINENO at command: $BASH_COMMAND && echo "exit status $exit_status" && exit' ERR
# if not 2 arguments passed, quit
if [ $# -ne 2 ]
then
    echo "arguments are missing to start script. 2 argumens are expected"
    exit 1
fi

# arguments
fastq_pass_dir_in="$1"
fastq_pass_dir_out="$2"

parent_directory=$(dirname "$fastq_pass_dir_out")

if [ ! -d "$parent_directory" ]; then
    echo "Parent directory for $fastq_pass_dir_out does not exist: $parent_directory"
    exit 1
fi

if [ ! -d "$fastq_pass_dir_out" ]; then
    mkdir "$fastq_pass_dir_out"
    echo "Directory created: $fastq_pass_dir_out"
else
    echo "Directory already exists: $fastq_pass_dir_out"
fi

for i in "$fastq_pass_dir_in"/barcode* ; do
    barcode_path=$(realpath "$i")
    echo "barcode path $barcode_path"
    barcode_name=$(basename "$barcode_path")
    echo "barcode_name $barcode_name"
    cat "$barcode_path/"*".fastq.gz" > "$fastq_pass_dir_out/${barcode_name}.fastq.gz"
done

