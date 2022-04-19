bigwigCompare=/media/stephen/Elements/envs/5353dfa1039e1a04b64654748c891041/bin/bigwigCompare

for p in 5 15 17 49 62 69 126; do
    for t in BEVA CTRL; do
        wiggletools median individual/PDOVCA_"${p}"_"${t}"-*.bigWig \
        | wiggletools bin 200 - \
        | grep -v GL000 | grep -v KI27 \
        | wigToBigWig -clip stdin ../genome/grch38.chrom.sizes PDOVCA_"${p}"_"${t}".bigWig
    done
    "${bigwigCompare}" \
        --bigwig1 PDOVCA_"${p}"_CTRL.bigWig --bigwig2 PDOVCA_"${p}"_BEVA.bigWig \
        --skipZeroOverZero --operation subtract --numberOfProcessors 4 \
        --outFileName PDOVCA_"${p}"_CTRL-minus-BEVA.bigWig
done

for t in BEVA CTRL; do

    wiggletools median \
        individual/PDOVCA_5_"${t}"-*.bigWig \
        individual/PDOVCA_126_"${t}"-*.bigWig \
        individual/PDOVCA_62_"${t}"-*.bigWig \
    | wiggletools bin 200 - \
    | grep -v GL000 | grep -v KI27 \
    | wigToBigWig -clip stdin ../genome/grch38.chrom.sizes responsive-"${t}".bigWig

    wiggletools median \
        individual/PDOVCA_15_"${t}"-*.bigWig \
        individual/PDOVCA_69_"${t}"-*.bigWig \
        individual/PDOVCA_17_"${t}"-*.bigWig \
        individual/PDOVCA_49_"${t}"-*.bigWig \
    | wiggletools bin 200 - \
    | grep -v GL000 | grep -v KI27 \
    | wigToBigWig -clip stdin ../genome/grch38.chrom.sizes resistant-"${t}".bigWig

    "${bigwigCompare}"  \
        --bigwig1 responsive-"${t}".bigWig --bigwig2 resistant-"${t}".bigWig \
        --skipZeroOverZero --operation subtract --numberOfProcessors 4 \
        --outFileName responsive-minus-resistant-"${t}".bigWig

done

for file in *-minus-*[AL].bigWig; do
     wiggletools "${file}" \
     | python scaleZ.py \
     | wigToBigWig -clip stdin ../genome/grch38.chrom.sizes "${file/.bigWig/-Z.bigWig}"
done
