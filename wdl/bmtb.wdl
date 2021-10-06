version 1.0

### bmtb.wdl ###
## Author: Charles Markello
## Description: Core black magic toolbox workflow for pedigree datasets

workflow bmtbWorkflow {
    meta {
        author: "Charles Markello"
        email: "cmarkell@ucsc.edu"
        description: "Core candidate variant analysis workflow for pedigree datasets."
    }
    input {
        # MISC inputs
        File INPUT_PED_FILE                         # Input trio .ped file for mosaicism detection
        File SNPEFF_DATABASE                        # Path to snpeff database .zip file for snpEff annotation functionality.
        
        # VCFtoShebang inputs
        File INPUT_VCF_FILE                         # Input cohort unrolled .vcf file
        String BYPASS = 'false'                     # Set to 'true' to bypass gender calling stage of vcftoshebang program
        Int CADD_LINES = 30000
        File CHROM_DIR
        File EDIT_DIR
        
        # CADD workflow inputs
        String GENOME_BUILD = "GRCh38"
        File CADD_DATA
        
        # BMTB workflow inputs
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
        Array[Int]+ AFFECTED_SIBLING_LIST           # Input list of sibling affected states (0 = affected, 1 = unaffected). Must follow same order as SAMPLE_NAME_SIBLING_LIST
    }
    
    # Run SNPEff Annotation
    call normalizeVCF as normalizeCohortVCF {
        input:
            in_sample_name=SAMPLE_NAME_SIBLING_LIST[0],
            in_bgzip_vcf_file=INPUT_VCF_FILE,
            in_small_resources=false
    }
    call snpEffAnnotateVCF as snpEffAnnotateCohortVCF {
        input:
            in_sample_name=SAMPLE_NAME_SIBLING_LIST[0],
            in_normalized_vcf_file=normalizeCohortVCF.output_normalized_vcf,
            in_snpeff_database=SNPEFF_DATABASE,
            in_small_resources=false
    }

    # Remove the decoy contig
    call run_remove_decoy_contigs {
        input:
            in_cohort_vcf=snpEffAnnotateCohortVCF.output_snpeff_annotated_vcf,
            in_genome_build=GENOME_BUILD
    }
    
    # Run mosaicism detection
    call run_detect_mosaicism {
        input:
            in_proband_name=SAMPLE_NAME_SIBLING_LIST[0],
            in_ped_file=INPUT_PED_FILE,
            in_cohort_vcf=snpEffAnnotateCohortVCF.output_snpeff_annotated_vcf
    }
    
    # Run VCFtoShebang conversion
    call run_vcf2shebang {
        input:
            in_proband_name=SAMPLE_NAME_SIBLING_LIST[0],
            in_cohort_vcf=run_remove_decoy_contigs.filtered_cohort_vcf,
            in_bypass=BYPASS,
            in_cadd_lines=CADD_LINES,
            in_chrom_file_dir=CHROM_DIR,
            in_edit_dir=EDIT_DIR
    }
    
    if(run_vcf2shebang.runCADD) {
        # Run CADD workflow
        call run_split_vcf {
            input:
                in_vcf=run_vcf2shebang.outputCADDVCF,
                in_split_lines=CADD_LINES
        }
        
        scatter (vcf_chunk_file in run_split_vcf.vcf_chunk_files) {
            call run_cadd {
                input:
                    in_chunk_vcf=vcf_chunk_file,
                    in_genome_build=GENOME_BUILD,
                    in_cadd_data=CADD_DATA
            }
        }
        
        Array[File?] cadd_annotated_vcf_chunks = run_cadd.outputCADD
        Array[File] cadd_annotated_vcf_chunks_list = select_all(cadd_annotated_vcf_chunks)
        
        call run_merge_annotated_vcf {
            input:
                in_cadd_output_chunks=cadd_annotated_vcf_chunks_list
        }
        
        call run_cadd_editor {
            input:
                in_vs_file=run_vcf2shebang.outputVS, 
                in_cadd_output=run_merge_annotated_vcf.merged_cadd_output_vcf
        }
    }
    
    # Run BMTB Candidate Analysis
    call run_bmtb {
        input:
            in_vs_file=select_first([run_cadd_editor.cadd_editor_output_vs, run_vcf2shebang.outputVS]),
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
        File output_mosaicism_report = run_detect_mosaicism.outputMOSAICISMREPORT
        File output_vs = run_bmtb.outputVS
    }
}

task normalizeVCF {
    input {
        String in_sample_name
        File in_bgzip_vcf_file
        Boolean in_small_resources
    }

    Int in_vgcall_cores = if in_small_resources then 6 else 6
    Int in_vgcall_disk = if in_small_resources then 1 else 25
    String in_vgcall_mem = if in_small_resources then "1" else "50"

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        bcftools norm -m-both --threads ~{in_vgcall_cores} -o ~{in_sample_name}.unrolled.vcf ~{in_bgzip_vcf_file}
    >>>
    output {
        File output_normalized_vcf = "~{in_sample_name}.unrolled.vcf"
    }
    runtime {
        preemptible: 2
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:95c212df20552fc74670d8f16d20099d9e76245eda6a1a6cfff4bd39e57be01b"
    }
}

task snpEffAnnotateVCF {
    input {
        String in_sample_name
        File in_normalized_vcf_file
        File? in_snpeff_database
        Boolean in_small_resources
    }

    Int in_vgcall_cores = if in_small_resources then 6 else 6
    Int in_vgcall_disk = if in_small_resources then 10 else 25
    String in_vgcall_mem = if in_small_resources then "10" else "50"

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        set -o xtrace
        #to turn off echo do 'set +o xtrace'

        unzip ~{in_snpeff_database}
        database_ref="GRCh38.99"
        if [[ "~{in_snpeff_database}" != *"GRCh38"* ]]; then
            database_ref="GRCh37.75"
        fi
        snpEff -Xmx40g -i VCF -o VCF -noLof -noHgvs -formatEff -classic -dataDir ${PWD}/data ${database_ref} ~{in_normalized_vcf_file} > ~{in_sample_name}.snpeff.unrolled.vcf
    >>>
    output {
        File output_snpeff_annotated_vcf = "~{in_sample_name}.snpeff.unrolled.vcf"
    }
    runtime {
        preemptible: 2
        cpu: in_vgcall_cores
        memory: in_vgcall_mem + " GB"
        disks: "local-disk " + in_vgcall_disk + " SSD"
        docker: "quay.io/biocontainers/snpeff:5.0--hdfd78af_1"
    }
}

task run_remove_decoy_contigs {
    input {
        File in_cohort_vcf
        String in_genome_build
    }
    String cohort_vcf_basename = basename(in_cohort_vcf)
    command <<<
        set -exu -o pipefail
        if [[ ~{in_genome_build} == *"GRCh38"* ]]; then
            vcftools --vcf ~{in_cohort_vcf} --not-chr hs38d1_decoys --recode-INFO-all --recode --out ~{cohort_vcf_basename}.filtered
        else
            vcftools --vcf ~{in_cohort_vcf} --not-chr hs37d5 --recode-INFO-all --recode --out ~{cohort_vcf_basename}.filtered
        fi
    >>>
    output {
        File filtered_cohort_vcf = "~{cohort_vcf_basename}.filtered.recode.vcf"
    }
    runtime {
        preemptible: 2
        memory: 20 + " GB"
        cpu: 1
        disks: "local-disk 20 SSD"
        docker: "biocontainers/vcftools:v0.1.16-1-deb_cv1"
    }
}

task run_detect_mosaicism {
    input {
        String in_proband_name
        File in_ped_file
        File in_cohort_vcf
    }
    String cohort_vcf_basename = basename(in_cohort_vcf)
    String pedfile_basename = basename(in_ped_file)
    command <<<
        set -exu -o pipefail
        
        cp /usr/src/app/phasing_config_file.txt .
        cp ~{in_ped_file} .
        cp ~{in_cohort_vcf} .
        sed -i "s|^PED_FILE.*|PED_FILE\t${PWD}/~{pedfile_basename}|" phasing_config_file.txt
        sed -i "s|^VCF_FILE.*|VCF_FILE\t${PWD}/~{cohort_vcf_basename}|" phasing_config_file.txt
        sed -i "s|^OUTPUT_FILE.*|OUTPUT_FILE\t${PWD}/~{in_proband_name}_mosaicism_output.txt|" phasing_config_file.txt
        sed -i "s|^PROBAND_NAME.*|PROBAND_NAME\t~{in_proband_name}|" phasing_config_file.txt
        python /usr/src/app/main.py -i phasing_config_file.txt
    >>>
    output {
        File outputMOSAICISMREPORT = "~{in_proband_name}_mosaicism_output.txt"
    }
    runtime {
        preemptible: 2
        memory: 100 + " GB"
        cpu: 1
        disks: "local-disk 100 SSD"
        docker: "quay.io/cmarkello/mosaicism_detector:latest"
    }
}

task run_vcf2shebang {
    input {
        String in_proband_name
        File in_cohort_vcf
        String in_bypass
        Int in_cadd_lines
        File in_chrom_file_dir
        File in_edit_dir
    }
    String cohort_vcf_basename = basename(in_cohort_vcf)
    String chrom_file_dir_basename = basename(in_chrom_file_dir, ".tar.gz")
    String edit_dir_basename = basename(in_edit_dir, ".tar.gz")
    command <<<
        set -exu -o pipefail
        
        tar -xzf ~{in_chrom_file_dir}
        tar -xzf ~{in_edit_dir}
        cp ~{in_cohort_vcf} .
        output_dir="vcf2shebang_output/"
        chrom_file_dir="$(basename ~{in_chrom_file_dir})"
        edit_dir="$(basename ~{in_edit_dir})"
        bypass_conf="NO"
        if [[ ~{in_bypass} == "true" ]]; then
            bypass_conf="YES"
        fi
        mkdir ${output_dir}
        cp /vcftoshebang/VCFtoShebang_Config.txt .
        sed -i "s|^PROBAND_NAME.*|PROBAND_NAME\t~{in_proband_name}|" VCFtoShebang_Config.txt
        sed -i "s|^OUTPUT_DIR.*|OUTPUT_DIR\t${output_dir}|" VCFtoShebang_Config.txt
        sed -i "s|^UNROLLED_VCF_PATH.*|UNROLLED_VCF_PATH\t~{cohort_vcf_basename}|" VCFtoShebang_Config.txt
        sed -i "s|^BYPASS.*|BYPASS\t${bypass_conf}|" VCFtoShebang_Config.txt
        sed -i "s|^CADD_LINES.*|CADD_LINES\t~{in_cadd_lines}|" VCFtoShebang_Config.txt
        sed -i "s|^CHROM_DIR.*|CHROM_DIR\t$PWD/~{chrom_file_dir_basename}|" VCFtoShebang_Config.txt
        sed -i "s|^EDIT_DIR.*|EDIT_DIR\t$PWD/~{edit_dir_basename}|" VCFtoShebang_Config.txt
        sed -i "s|^EDITOR_CONFIG.*|EDITOR_CONFIG\t/vcftoshebang/edit_config.txt|" VCFtoShebang_Config.txt

        java -XX:+UnlockExperimentalVMOptions -XX:ActiveProcessorCount=8 -cp /vcftoshebang/VCFtoShebang.jar:/vcftoshebang/json_simple.jar Runner VCFtoShebang_Config.txt
        
        # Check for CADD input
        echo "false" > vcf2shebang_output/cadd_input_available
        while IFS= read -r line
        do
            if [[ line != *"#"* ]]; then
                echo "true" > vcf2shebang_output/cadd_input_available
            fi
        done < "vcf2shebang_output/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_Input_Files/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt.gz"
    >>>
    output {
        File outputVS = "vcf2shebang_output/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_shebang.vs"
        File outputCADDVCF = "vcf2shebang_output/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_Input_Files/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt.gz"
        Boolean runCADD = read_boolean("vcf2shebang_output/cadd_input_available")
    }
    runtime {
        preemptible: 2
        memory: 100 + " GB"
        cpu: 8
        disks: "local-disk 400 SSD"
        docker: "quay.io/cmarkello/vcf2shebang_grch38:latest"
    }
}

task run_split_vcf {
    input {
        File in_vcf
        Int in_split_lines
    }
    String vcf_basename = basename(in_vcf, ".txt.gz")
    command <<<
        set -exu -o pipefail
        
        tempname="~{vcf_basename}_tmp"
        zcat ~{in_vcf} | grep -v '^GL|^#' | cut -f 1-5 | split -l ~{in_split_lines} --additional-suffix=".vcf" -d - ${tempname}
        ls -l
    >>>
    output {
        Array[File] vcf_chunk_files = glob(".*tmp.*.vcf")
    }
    runtime {
        preemptible: 2
        memory: 20 + " GB"
        cpu: 4
        disks: "local-disk 20 SSD"
        docker: "ubuntu@sha256:2695d3e10e69cc500a16eae6d6629c803c43ab075fa5ce60813a0fc49c47e859"
    }
}

task run_cadd {
    input {
        File in_chunk_vcf
        String in_genome_build
        File in_cadd_data
    }
    String cadd_data_dir_basename = basename(in_cadd_data, ".tar.gz")
    command <<<
        set -exu -o pipefail
        tar -xzf ~{in_cadd_data}
        base=$(basename -s .vcf "~{in_chunk_vcf}")
        /bin/bash /usr/src/app/CADD.sh -v "v1.6" -g ~{in_genome_build} -o ${base}_out.tsv.gz -d ~{cadd_data_dir_basename} ~{in_chunk_vcf}
    >>>
    output {
        File outputCADD = glob("*_out.tsv.gz")[0]
    }
    runtime {
        preemptible: 2
        memory: 100 + " GB"
        cpu: 6
        disks: "local-disk 150 SSD"
        docker: "quay.io/cmarkello/cadd_1.6:latest"
    }
}

task run_merge_annotated_vcf {
    input {
        Array[File] in_cadd_output_chunks
    }
    command <<<
        set -exu -o pipefail
        
        
        for cadd_output_chunk in $(ls ~{sep=" " in_cadd_output_chunks} | sort) ; do
            zcat ${cadd_output_chunk} >> merged_CADDv1.6_offline_unsorted
        done
        sort -k1,1 -k2,2n merged_CADDv1.6_offline_unsorted > merged_CADDv1.6_offline.vcf && \
        rm -f merged_CADDv1.6_offline_unsorted && \
        python3 /usr/src/app/CADD_offline_mito_postprocessing.py -c /usr/src/app/whole_mito_SNP_pp2_predictions_sorted.txt -i merged_CADDv1.6_offline.vcf -o merged_CADDv1.6_offline_proper_format.vcf && \
        rm -f merged_CADDv1.6_offline.vcf
    >>>
    output {
        File merged_cadd_output_vcf = "merged_CADDv1.6_offline_proper_format.vcf"
    }
    runtime {
        preemptible: 2
        memory: 50 + " GB"
        cpu: 4
        disks: "local-disk 100 SSD"
        docker: "quay.io/cmarkello/cadd_1.6:latest"
    }
}

task run_cadd_editor {
    input {
        File in_vs_file
        File in_cadd_output
    }
    command <<<
        set -exu -o pipefail
        
        java -cp /cadd_edit/NewCaddEditor.jar:/cadd_edit/commons-cli-1.4.jar NewCaddEditor --input_vs ~{in_vs_file} --output_cadd ~{in_cadd_output} --output_vs "cadd_editor_output.vs"
    >>>
    output {
        File cadd_editor_output_vs = "cadd_editor_output.vs"
    }
    runtime {
        preemptible: 2
        memory: 20 + " GB"
        cpu: 2
        disks: "local-disk 100 SSD"
        docker: "quay.io/cmarkello/cadd_editor:latest"
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
        
        rm ${PWD}/Configs/BAM_Directory_Config.txt
        touch ${PWD}/Configs/BAM_Directory_Config.txt
        echo -e "~{in_maternal_id}\t${PWD}/input_maternal.bam" >> ${PWD}/Configs/BAM_Directory_Config.txt
        echo -e "~{in_paternal_id}\t${PWD}/input_paternal.bam" >> ${PWD}/Configs/BAM_Directory_Config.txt 
         
        sibling_id_list=(~{sep=" " in_sibling_ids})
        sibling_gender_list=(~{sep=" " in_sibling_genders})
        sibling_affected_list=(~{sep=" " in_sibling_affected})
        sibling_bam_list=(~{sep=" " in_sibling_bams})
        sibling_bam_index_list=(~{sep=" " in_sibling_bams_index})
        
        SIB_IDS_STRING=""
        
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
                sed -i "s|^PB_ID.*|PB_ID\t${SIBLING_ID}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
                sed -i "s|^PB_GENDER.*|PB_GENDER\t${SIBLING_GENDER}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
                if [[ SIBLING_GENDER == *"0"* ]]; then
                    sed -i "s|^FILTER_HEMI.*|FILTER_HEMI\t1|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
                else
                    sed -i "s|^FILTER_HEMI.*|FILTER_HEMI\t0|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
                fi
            else
                SIB_IDS_STRING="${SIB_IDS_STRING}${SIBLING_ID},${SIBLING_GENDER},${SIBLING_AFFECTED};"
            fi
        done
        
        sed -i "s|^OUT_FILE_NAME.*|OUT_FILE_NAME\t~{proband_id}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|^VS_FILE_PATH.*|VS_FILE_PATH\t${PWD}/input_vs.vs|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|^CONFIG_FILE_DIREC.*|CONFIG_FILE_DIREC\t${PWD}/Configs|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|^FATHER_ID.*|FATHER_ID\t~{in_paternal_id}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|^MOTHER_ID.*|MOTHER_ID\t~{in_maternal_id}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        sed -i "s|^SIB_IDS.*|SIB_IDS\t${SIB_IDS_STRING%?}|" ${PWD}/Configs/BMTB_Genome_Input_Config.txt
         
        java -cp /bmtb/bmtb.jar:/bmtb/htsjdk-2.19.0-47-gc5ed6b7-SNAPSHOT.jar general.Runner ${PWD}/Configs/BMTB_Genome_Input_Config.txt
        tar -czvf "~{proband_id}_BlackBox_Output.tar.gz" "~{proband_id}_BlackBox_Output"
    >>>
    output {
        File outputVS = "~{proband_id}_BlackBox_Output.tar.gz"
    }
    runtime {
        preemptible: 2
        memory: 100 + " GB"
        cpu: 6
        disks: "local-disk 150 SSD"
        docker: "quay.io/cmarkello/bmtb_grch38:latest"
    }
}

