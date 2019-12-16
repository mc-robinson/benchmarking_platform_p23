# An example function call used to create the rankings
# python calculate_scored_lists.py -n 10 -f ./fp_file.txt -s Tanimoto -o ./test_results/
# should be run in p2-env

# for LR
python calculate_scored_lists_LR.py -n 10 -f rdk6 -s Tanimoto -o ./test_results/

# Note in the RF code, I had to disable his special parallel build function

# Note that ap fp is failing with RF, just `killed` message nothing else
