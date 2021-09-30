version development

### bmtb.wdl ###
## Author: Charles Markello
## Description: Core black magic toolbox workflow for pedigree datasets

workflow bmtbWorkflow {
    input {
        # BMTB workflow inputs
        File INPUT_VS_FILE                          # Input cohort processed varsifter file
        File MATERNAL_INPUT_BAM_FILE                # Input maternal .bam file
        File MATERNAL_INPUT_BAM_FILE_INDEX          # Input maternal .bai index of .bam file
        File PATERNAL_INPUT_BAM_FILE                # Input paternal .bam file
        File PATERNAL_INPUT_BAM_FILE_INDEX          # Input paternal .bai index of .bam file
        Array[File]+ SIBLING_BAM_FILE_LIST          # Input list of sibling .bam files. Proband must be first in list.
        Array[File]+ SIBLING_BAM_FILE_INDEX_LIST    # Input list of .bai indices of .bam files. Must follow same sample order as SIBLING_BAM_FILE_LIST.
        String SAMPLE_NAME_MATERNAL                 # Sample name for the mother
        String SAMPLE_NAME_PATERNAL                 # Sample name for the father
        Array[String]+ SAMPLE_NAME_SIBLING_LIST     # Input list of sibling sample ids. Must follow same order as SIBLING_BAM_FILE_LIST.
        Array[Int]+ GENDER_SIBLING_LIST             # Input list of sibling gender states (0 = male, 1 = female). Must follow same order as SAMPLE_NAME_SIBLING_LIST
        Array[Int]+ AFFECTED_SIBLING_LIST           # Input list of sibling affected states (0 = unaffected, 1 = affected). Must follow same order as SAMPLE_NAME_SIBLING_LIST
    }
    
    # Run BMTB Candidate Analysis
    call run_bmtb {
        input:
            in_vs_file=INPUT_VS_FILE,
            in_maternal_bam=MATERNAL_INPUT_BAM_FILE,
            in_maternal_bam_index=MATERNAL_INPUT_BAM_FILE_INDEX,
            in_paternal_bam=PATERNAL_INPUT_BAM_FILE,
            in_paternal_bam_index=PATERNAL_INPUT_BAM_FILE_INDEX,
            in_sibling_bams=SIBLING_BAM_FILE_LIST,
            in_sibling_bams_index=SIBLING_BAM_FILE_INDEX_LIST,
            in_maternal_id=SAMPLE_NAME_MATERNAL,
            in_paternal_id=SAMPLE_NAME_PATERNAL,
            in_sibling_ids=SAMPLE_NAME_SIBLING_LIST,
            in_sibling_genders=GENDER_SIBLING_LIST,
            in_sibling_affected=AFFECTED_SIBLING_LIST
    }
    
    output {
        File output_vs = run_bmtb.outputVS
    }
}

task run_bmtb {
    input {
        File in_vs_file
        File in_maternal_bam
        File in_maternal_bam_index
        File in_paternal_bam
        File in_paternal_bam_index
        Array[File]+ in_sibling_bams
        Array[File]+ in_sibling_bams_index
        String in_maternal_id
        String in_paternal_id
        Array[String]+ in_sibling_ids
        Array[Int]+ in_sibling_genders
        Array[Int]+ in_sibling_affected
    }
    
    Int num_siblings = length(in_sibling_ids)
    String proband_id = in_sibling_ids[0]
    
    command <<<
        set -exu -o pipefail
        
        ln -s ~{in_maternal_bam} input_maternal.bam
        ln -s ~{in_maternal_bam_index} input_maternal.bam.bai
        ln -s ~{in_paternal_bam} input_paternal.bam
        ln -s ~{in_paternal_bam_index} input_paternal.bam.bai
        ln -s ~{in_vs_file} input_vs.vs
        
        cp -r /bmtb/Configs ${PWD}/Configs
        rm -f ${PWD}/Configs/BAM_Directory_Config.txt
        touch ${PWD}/Configs/BAM_Directory_Config.txt
        echo -e "~{in_maternal_id}\t${PWD}/input_maternal.bam" >> ${PWD}/Configs/BAM_Directory_Config.txt
        echo -e "~{in_paternal_id}\t${PWD}/input_paternal.bam" >> ${PWD}/Configs/BAM_Directory_Config.txt 
         
        sibling_id_list=(~{sep=" " in_sibling_ids})
        sibling_gender_list=(~{sep=" " in_sibling_genders})
        sibling_affected_list=(~{sep=" " in_sibling_affected})
        sibling_bam_list=(~{sep=" " in_sibling_bams})
        sibling_bam_index_list=(~{sep=" " in_sibling_bams_index})
        
        SIB_IDS_STRING="NA"
        
        for (( i=0; i<~{num_siblings}; i++ )) ; do
            SIBLING_ID=${sibling_id_list[$i]}
            SIBLING_GENDER=${sibling_gender_list[$i]}
            SIBLING_AFFECTED=${sibling_affected_list[$i]}
            SIBLING_BAM_FILE=${sibling_bam_list[$i]}
            SIBLING_BAM_INDEX_FILE=${sibling_bam_index_list[$i]}
            ln -s ${SIBLING_BAM_FILE} input_${SIBLING_ID}.bam
            ln -s ${SIBLING_BAM_INDEX_FILE} input_${SIBLING_ID}.bam.bai
            echo -e "${SIBLING_ID}\t${PWD}/input_${SIBLING_ID}.bam" >> ${PWD}/Configs/BAM_Directory_Config.txt
            if [ ${i} == 0 ]; then
                sed -i "s|.*PB_ID.*|PB_ID\t${SIBLING_ID}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
                sed -i "s|.*PB_GENDER.*|PB_GENDER\t${SIBLING_GENDER}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
            else
                if [ ${SIB_IDS_STRING} == "NA" ]; then
                    SIB_IDS_STRING="${SIBLING_ID},${SIBLING_GENDER},${SIBLING_AFFECTED};"
                else
                    SIB_IDS_STRING="${SIB_IDS_STRING}${SIBLING_ID},${SIBLING_GENDER},${SIBLING_AFFECTED};"
                fi
            fi
        done
        
        sed -i "s|.*OUT_FILE_NAME.*|OUT_FILE_NAME\t~{proband_id}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|.*VS_FILE_PATH.*|VS_FILE_PATH\t${PWD}/input_vs.vs|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|.*CONFIG_FILE_DIREC.*|CONFIG_FILE_DIREC\t${PWD}/Configs|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|.*FATHER_ID.*|FATHER_ID\t~{in_paternal_id}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|.*MOTHER_ID.*|MOTHER_ID\t~{in_maternal_id}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|.*SIB_IDS.*|SIB_IDS\t${SIB_IDS_STRING%?}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
         
        java -cp /bmtb/bmtb.jar:/bmtb/htsjdk-2.19.0-47-gc5ed6b7-SNAPSHOT.jar general.Runner ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        tar -czvf "~{proband_id}_BlackBox_Output.tar.gz" "~{proband_id}_BlackBox_Output"
    >>>
    output {
        File outputVS = "~{proband_id}_BlackBox_Output.tar.gz"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 6
        docker: "quay.io/cmarkello/bmtb:latest"
    }
}

