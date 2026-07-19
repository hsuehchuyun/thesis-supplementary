#!/bin/bash 
#SBATCH -A MST114487 
#SBATCH -J ldsc_EUR
#SBATCH -c 4 
#SBATCH -p ct56
#SBATCH --mem=32G
#SBATCH -o eur_out.log 
#SBATCH -e eur_err.log 

module load old-module
ml load pkg/Anaconda3
source /opt/ohpc/Taiwania3/pkg/anaconda3.4.9.2/etc/profile.d/conda.sh
conda activate ldsc
# my_env

dir0=/work/u7080475/0_program/1000_genome/eur_w_ld_chr
dir1=/work/u7080475/0_formatdata
dir3=/work/u7080475/0_program/1000_genome
ldsc_py=/work/u7080475/0_program/ldsc


# "phenotype_name:sample_size"
declare -a PHENOS=(
    "8_mdd_howD19:500199"
    "7_sczeur_vasT22:130644"
    "TRD_CAD_ukbcard:640653"
    "40_t2d_mahA18:898130"
    "TRD_crp_sarS22:575531"
    "TRD_DrnkWkEur_greS22:666978"
    "TRD_AgeSmkEur_greS22:323386"
    "TRD_CigDayEur_greS22:326497"
    "TRD_SmkInitEur_greS22:805431"
    "TRD_SmkCesEur_gre22:388313"
    "TRD_EA_jamL18:766345"
    "Inflam_neut_minC20:519288"
    "Inflam_lymp_minC20:524923"
    "TRD_strokeeur_aniM22:1308460"
    "TRD_cardstrokeeur_aniM22:1245612"
    "TRD_iscstrokeeur_aniM22:1296908"
    "TRD_svesstrokeeur_aniM22:1241619"
    "TRD_lvesstrokeeur_aniM22:1241207"
    "TRD_HtEur_loiY22:1597374"
    "TRD_wt_ukbv3:360116"
    "37_bmi_sarP19:806834"
    "TRD_bfr_ukbv3:354628"
    "TRD_wc_ukbv3:360564"
    "TRD_hc_ukbv3:360521"
    "38_whr_sarP19:697734"
    "TRD_dbp_ukbv3:340162"
    "TRD_sbp_ukbv3:340159"
    "TRD_HrtRate_zhaZ19:458969"
    "TRD_wbc_ukbv3:350470"
    "TRD_rbc_ukbv3:350475"
    "TRD_hb_ukbv3:350474"
    "TRD_hct_ukbv3:350475"
    "TRD_plt_ukbv3:350474"
    "TRD_bun_ukbv3:344052"
    "TRD_cr_ukbv3:344104"
    "TRD_malb_ukbv3:108706"
    "TRD_ua_ukbv3:343836"
    "TRD_tbil_ukbv3:342829"
    "TRD_alt_ukbv3:344136"
    "TRD_ast_ukbv3:342990"
    "TRD_ggt_ukbv3:344104"
    "TRD_alb_ukbv3:315268"
    "TRD_fev_ukbv3:329404"
    "TRD_fvc_ukbv3:329404"
    "TRD_fg_ukbv3:314916"
    "TRD_hba1c_ukbv3:344182"
    "TRD_tc_ukbv3:344278"
    "TRD_hdl_ukbv3:315133"
    "TRD_ldl_ukbv3:343621"
    "TRD_tg_ukbv3:343992"
)

cd /home/u7080475/ldsc

mkdir -p munge/EUR
mkdir -p h2/EUR
mkdir -p gc/EUR

# Munge 
echo "==================== Starting Munge ===================="
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
        --out munge/EUR/${pheno_name} \
        --chunksize 500000 \
        --ignore variant,minor_allele,maf,low_confidence_variant,info,SNPID,chrpos,HWEP,AC,ytx,tstat
done

# h2
echo "==================== Starting h2 Calculation ===================="
for pheno_info in "${PHENOS[@]}"; do
    pheno_name="${pheno_info%%:*}"
    
    echo "Calculating h2 for ${pheno_name}..."
    python -u $ldsc_py/ldsc.py \
        --h2 munge/EUR/${pheno_name}.sumstats.gz \
        --ref-ld-chr ${dir0}/ \
        --w-ld-chr ${dir0}/ \
        --out h2/EUR/h2_${pheno_name}
done

# rg
echo "==================== Starting rg Calculation ===================="
total_phenos=${#PHENOS[@]}
for ((i=0; i<${total_phenos}; i++)); do
    pheno1_name="${PHENOS[$i]%%:*}"
    
    for ((j=i+1; j<${total_phenos}; j++)); do
        pheno2_name="${PHENOS[$j]%%:*}"
        
        echo "Calculating rg between ${pheno1_name} and ${pheno2_name}..."
        python -u $ldsc_py/ldsc.py \
            --rg munge/EUR/${pheno1_name}.sumstats.gz,munge/EUR/${pheno2_name}.sumstats.gz \
            --ref-ld-chr ${dir0}/ \
            --w-ld-chr ${dir0}/ \
            --out gc/EUR/rg_EUR_${pheno1_name}_${pheno2_name}
    done
done

echo "==================== EUR Task finished at $(date) ===================="



