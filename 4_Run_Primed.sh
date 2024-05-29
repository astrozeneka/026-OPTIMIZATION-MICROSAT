#!/bin/bash

if [ $# -ne 3 ]; then
	echo "Usage: $0 <marker_number> <replicate_id> <pheromone npy file>"
	exit 1
fi

source ~/.bashrc
conda activate Rscript
R_LIBS=/home/kornsorn02/ryan-agb/Projects/ACO_MicrosatelliteMarkers/r_lib

start_time=$(date +%s)
echo "Script started at: $(date)"
python3 4_Optimize.py -f data/gallus-gallus.str -i 652 -l 28 -m "$1" -s "$2" -e 120 -p "$3" \
	-b /home/kornsorn02/.conda/envs/Rscript/bin/Rscript
end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Script finished at: $(date)"
echo "Total execution time: $execution_time seconds"


