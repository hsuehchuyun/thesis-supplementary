#!/bin/bash 
#SBATCH -A MST114487 
#SBATCH -J ldsc_EAS
#SBATCH -c 4 
#SBATCH -p ct56
#SBATCH --mem=32G
#SBATCH -o eas_out.log 
#SBATCH -e eas_err.log 

module load old-module
ml load pkg/Anaconda3
source /opt/ohpc/Taiwania3/pkg/anaconda3.4.9.2/etc/profile.d/conda.sh
conda activate ldsc

dir0=/work/u7080475/0_program/1000_genome/eas_ldscores
dir1=/work/u7080475/0_formatdata
dir3=/work/u7080475/0_program/1000_genome
ldsc_py=/work/u7080475/0_program/ldsc

# "phenotype_name:sample_size"
declare -a PHENOS=(
    "TRD_mddasi_olgG21:169805"
    "TRD_sczasn_maxL19:58140"
    "TRD_cadasn_satK20:168228"
    "TRD_t2d_twbasn:538294"
    "TRD_crpasn_saoS21:83025"
    "TRD_DrnkWkEas_greS22:90852"
    "TRD_AgeSmkEas_greS22:59148"
    "TRD_CigDayEas_greS22:104664"
    "TRD_SmkInitEas_greS22:217901"
    "TRD_SmkCesEas_gre22:106992"
    "TRD_EA_twb:104722"
    "Inflam_neutasia_minC20:78744"
    "Inflam_lympasia_minC20:89266"
    "TRD_strokeeas_aniM22:264655"
    "TRD_cardstrokeeas_aniM22:238168"
    "TRD_iscstrokeeas_aniM22:256274"
    "TRD_svesstrokeeas_aniM22:242774"
    "TRD_lvesstrokeeas_aniM22:238977"
    "TRD_ht_twbasn:456471"
    "TRD_wt_twbbbj:258034"
    "TRD_bmi_twbbbj:256450"
    "TRD_bfr_twb:92615"
    "TRD_wc_twb:92615"
    "TRD_hc_twb:92615"
    "TRD_whr_twb:92615"
    "TRD_dbp_twbbbj:238130"
    "TRD_sbp_twbbbj:238120"
    "TRD_hr_twb:92615"
    "TRD_wbc_twbbbj:246970"
    "TRD_rbc_twbbbj:246127"
    "TRD_hb_twbbbj:245062"
    "TRD_hct_twbbbj:245630"
    "TRD_plt_twbbbj:241238"
    "TRD_bun_twbbbj:241382"
    "TRD_cr_twbbbj:242881"
    "TRD_malb_twb:92615"
    "TRD_ua_twbbbj:222020"
    "TRD_tbil_twbbbj:216956"
    "TRD_alt_twbbbj:243160"
    "TRD_ast_twbbbj:242683"
    "TRD_ggt_twbbbj:226086"
    "TRD_alb_twbbbj:213154"
    "TRD_fev_twb:62901"
    "TRD_fvc_twb:62901"
    "TRD_fg_twbbbj:225951"
    "TRD_hba1c_twbbbj:163836"
    "TRD_tc_twbbbj:228423"
    "TRD_hdl_twbbbj:167585"
    "TRD_ldl_twbbbj:165481"
    "TRD_tg_twbbbj:204282"
)

cd /home/u7080475/ldsc

mkdir -p munge/EAS
mkdir -p h2/EAS
mkdir -p gc/EAS

# Munge
echo "==================== Starting Munge (EAS) ===================="
for pheno_info in "${PHENOS[@]}"; do
    pheno_name="${pheno_info%%:*}"
    sample_size="${pheno_info##*:}"
    
    echo "Munging ${pheno_name} (N=${sample_size})..."
    python -u $ldsc_py/munge_sumstats.py \
        --sumstats ${dir1}/${pheno_name}.gz \
        --snp SNP \
        --a1 effect_allele \
        --a2 other_allele \
        --frq eaf \
        --p pval \
        --N ${sample_size} \
        --signed-sumstats beta,0 \
        --merge-alleles $dir3/w_hm3.snplist \
        --out munge/EAS/${pheno_name} \
        --chunksize 500000
done

# h2
echo "==================== Starting h2 Calculation (EAS) ===================="
for pheno_info in "${PHENOS[@]}"; do
    pheno_name="${pheno_info%%:*}"
    
    echo "Calculating h2 for ${pheno_name}..."
    python -u $ldsc_py/ldsc.py \
        --h2 munge/EAS/${pheno_name}.sumstats.gz \
        --ref-ld-chr ${dir0}/ \
        --w-ld-chr ${dir0}/ \
        --out h2/EAS/h2_${pheno_name}
done

# rg
echo "==================== Starting rg Calculation (EAS) ===================="
total_phenos=${#PHENOS[@]}
for ((i=0; i<${total_phenos}; i++)); do
    pheno1_name="${PHENOS[$i]%%:*}"
    
    for ((j=i+1; j<${total_phenos}; j++)); do
        pheno2_name="${PHENOS[$j]%%:*}"
        
        echo "Calculating rg between ${pheno1_name} and ${pheno2_name}..."
        python -u $ldsc_py/ldsc.py \
            --rg munge/EAS/${pheno1_name}.sumstats.gz,munge/EAS/${pheno2_name}.sumstats.gz \
            --ref-ld-chr ${dir0}/ \
            --w-ld-chr ${dir0}/ \
            --out gc/EAS/rg_EAS_${pheno1_name}_${pheno2_name}
    done
done

echo "==================== EAS Task finished at $(date) ===================="