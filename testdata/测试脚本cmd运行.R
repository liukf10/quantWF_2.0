. ~/conda.sh #进入conda
source activate R4.0 

cd ~/testdata #go to the working directory

Rscript /mnt/sharing/quantWF_2.0/step2_batch_adj_R \
--quantfile "test_for_batch_adj.RData" --output "batchtest" --funcdir "/mnt/sharing/quantWF_2.0/" \
--keepvar "Group,Sex" --batchname "Study" --evaluate_method "pca&pvca&mds&bic"

#help
Rscript /mnt/sharing/quantWF_2.0/step0_MS_R --help

Rscript /mnt/sharing/quantWF_2.0/step0_MS_R \
--quantfile "phos_data.txt" --output "phos" --IDinfCol "1,2,3,4,5,6,7"

