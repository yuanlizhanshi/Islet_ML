macs2 callpeak -t SRR9336443_GSM3900796_Human_islet_control_rep_HI22_ChIP-Seq_Homo_sapiens_ChIP-Seq_rmdup.bam \
-f BAM -g hs -n H3K27ac_HI22_islet -q 0.01 --outdir ../peaks
macs2 callpeak -t SRR9336444_GSM3900797_Human_islet_control_rep_HI19_ChIP-Seq_Homo_sapiens_ChIP-Seq_rmdup.bam \
-f BAM -g hs -n H3K27ac_HI19_islet -q 0.01 --outdir ../peaks
macs2 callpeak -t SRR9336445_GSM3900798_Human_islet_control_rep_HI37_ChIP-Seq_Homo_sapiens_ChIP-Seq_rmdup.bam \
-f BAM -g hs -n H3K27ac_HI37_islet -q 0.01 --outdir ../peaks
macs2 callpeak -t SRR9336446_GSM3900799_Human_islet_control_rep_HI40_ChIP-Seq_Homo_sapiens_ChIP-Seq_rmdup.bam \
-f BAM -g hs -n H3K27ac_HI40_islet -q 0.01 --outdir ../peaks
