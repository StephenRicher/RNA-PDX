for p in 5 15 17 49 62 69 126; do
    for t in BEVA CTRL; do
        wiggletools mean individual/PDOVCA_"${p}"_"${t}"-*.bigWig \
        | wiggletools bin 100 - \
        | grep -v GL000 | grep -v KI27 \
        | wigToBigWig -clip stdin ../genome/grch38.chrom.sizes PDOVCA_"${p}"_"${t}".bigWig
    done
done

for t in BEVA CTRL; do
    wiggletools mean individual/PDOVCA_*_"${t}"-*.bigWig \
    | wiggletools bin 100 - \
    | grep -v GL000 | grep -v KI27 \
    | wigToBigWig -clip stdin ../genome/grch38.chrom.sizes "${t}".bigWig
done

wiggletools mean \
    individual/PDOVCA_5_*-*.bigWig \
    individual/PDOVCA_126_*-*.bigWig \
    individual/PDOVCA_62_*-*.bigWig \
| wiggletools bin 100 - \
| grep -v GL000 | grep -v KI27 \
| wigToBigWig -clip stdin ../genome/grch38.chrom.sizes responsive.bigWig

wiggletools mean \
    individual/PDOVCA_15_*-*.bigWig \
    individual/PDOVCA_69_*-*.bigWig \
    individual/PDOVCA_17_*-*.bigWig \
    individual/PDOVCA_49_*-*.bigWig \
| wiggletools bin 100 - \
| grep -v GL000 | grep -v KI27 \
| wigToBigWig -clip stdin ../genome/grch38.chrom.sizes resistant.bigWig
