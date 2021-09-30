version 1.0

### vcftoshebang.wdl ###
## Author: Charles Markello
## Description: Core black magic toolbox workflow for pedigree datasets

workflow bmtbWorkflow {
    input {
        # VCFtoShebang inputs
        File INPUT_VCF_FILE                         # Input cohort unrolled .vcf file
        String SAMPLE_NAME
        String BYPASS
        Int CADD_LINES
        String CHROM_DIR
        String EDIT_DIR
        
        # CADD workflow inputs
        Int SPLIT_LINES
        String GENOME_BUILD
        String CADD_DATA
    }
    
    # Run VCFtoShebang conversion
    call run_vcf2shebang {
        input:
            in_proband_name=SAMPLE_NAME,
            in_cohort_vcf=INPUT_VCF_FILE,
            in_bypass=BYPASS,
            in_cadd_lines=CADD_LINES,
            in_chrom_file_dir=CHROM_DIR,
            in_edit_dir=EDIT_DIR
    }
   
    ## TODO
    if(length(run_vcf2shebang.outputCADDVCF) > 0) {
        # Run CADD workflow
        call run_split_vcf {
            input:
                in_vcf=run_vcf2shebang.outputCADDVCF[0],
                in_split_lines=SPLIT_LINES
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
                in_vs_file=run_vcf2shebang.outputVS, # Not sure where this file comes from
                in_cadd_output=run_merge_annotated_vcf.merged_cadd_output_vcf
        }
    }
    output {
        File output_vs = select_first([run_cadd_editor.cadd_editor_output_vs, run_vcf2shebang.outputVS])
    }
}

task run_vcf2shebang {
    input {
        String in_proband_name
        File in_cohort_vcf
        String in_bypass
        Int in_cadd_lines
        String in_chrom_file_dir
        String in_edit_dir
    }
    command <<<
        set -exu -o pipefail
        
        output_dir="vcf2shebang_output/"
        chrom_file_dir="$(basename ~{in_chrom_file_dir})"
        edit_dir="$(basename ~{in_edit_dir})"
        glchrom_config_filename="glchrom_config.txt"
        mkdir ${output_dir}
        cp /vcftoshebang/VCFtoShebang_Config.txt .
        sed -i "s|.*PROBAND_NAME.*|PROBAND_NAME\t~{in_proband_name}|" VCFtoShebang_Config.txt
        sed -i "s|.*OUTPUT_DIR.*|OUTPUT_DIR\t${output_dir}|" VCFtoShebang_Config.txt
        sed -i "s|.*UNROLLED_VCF_PATH.*|UNROLLED_VCF_PATH\t~{in_cohort_vcf}|" VCFtoShebang_Config.txt
        sed -i "s|.*BYPASS.*|BYPASS\t~{in_bypass}|" VCFtoShebang_Config.txt
        sed -i "s|.*CADD_LINES.*|CADD_LINES\t~{in_cadd_lines}|" VCFtoShebang_Config.txt
        sed -i "s|.*CHROM_DIR.*|CHROM_DIR\t$PWD/${chrom_file_dir}|" VCFtoShebang_Config.txt
        sed -i "s|.*EDIT_DIR.*|EDIT_DIR\t$PWD/${edit_dir}|" VCFtoShebang_Config.txt
        ln -s /vcftoshebang/glchrom_config.txt .
        sed -i "s|.*GL_CONFIG.*|GL_CONFIG\t$PWD/${glchrom_config_filename}|" VCFtoShebang_Config.txt
        
        java -XX:+UnlockExperimentalVMOptions -XX:ActiveProcessorCount=32 -cp /vcftoshebang/VCFtoShebang.jar:/vcftoshebang/json_simple.jar Runner VCFtoShebang_Config.txt && \
        if [ $(cat vcf2shebang_output/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_Input_Files/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt | grep -v "^#" | wc -l) -gt 0 ] ; then
            mv vcf2shebang_output/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_Input_Files/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt.gz ~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt.gz;
        fi
    >>>
    output {
        File outputVS = "vcf2shebang_output/~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_shebang.vs"
        Array[File] outputCADDVCF = glob("~{in_proband_name}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt.gz")
    }
    runtime {
        memory: 100 + " GB"
        cpu: 32
        docker: "quay.io/cmarkello/vcf2shebang:latest"
        mount1: in_chrom_file_dir
        mount2: in_edit_dir
    }
}

task run_split_vcf {
    input {
        File in_vcf
        Int in_split_lines
    }
    command <<<
        set -exu -o pipefail
        
        base=$(basename -s .txt.gz "~{in_vcf}")
        tempname="${base}_tmp"
        zcat ~{in_vcf} | grep -v "^GL\|^#" | cut -d$'\t' -f 1-5 | split -l ~{in_split_lines} --additional-suffix=".vcf" -d - ${tempname}
    >>>
    output {
        Array[File] vcf_chunk_files = glob("*tmp*.vcf")
    }
    runtime {
        memory: 20 + " GB"
        cpu: 4
        docker: "ubuntu@sha256:2695d3e10e69cc500a16eae6d6629c803c43ab075fa5ce60813a0fc49c47e859"
    }
}

task run_cadd {
    input {
        File in_chunk_vcf
        String in_genome_build
        String in_cadd_data
    }
    command <<<
        set -exu -o pipefail
        
        source activate $(head -1 /usr/src/app/environment.yml | cut -d' ' -f2)
        base=$(basename -s .vcf "~{in_chunk_vcf}")
        cadd_data_dir_base=$(basename "~{in_cadd_data}")
        /bin/bash /usr/src/app/CADD.sh -v "v1.5" -g ~{in_genome_build} -o $PWD/${base}_out.tsv.gz -d $PWD/${cadd_data_dir_base} ~{in_chunk_vcf}
    >>>
    output {
        File outputCADD = glob("*_out.tsv.gz")[0]
    }
    runtime {
        memory: 100 + " GB"
        cpu: 6
        docker: "quay.io/cmarkello/cadd:latest"
        mount1: in_cadd_data
    }
}

task run_merge_annotated_vcf {
    input {
        Array[File] in_cadd_output_chunks
    }
    command <<<
        set -exu -o pipefail
        
        source activate $(head -1 /usr/src/app/environment.yml | cut -d' ' -f2)
        for cadd_output_chunk in $(ls ~{sep=" " in_cadd_output_chunks} | sort) ; do
            zcat ${cadd_output_chunk} >> merged_CADDv1.5_offline_unsorted
        done
        sort -k1,1 -k2,2n merged_CADDv1.5_offline_unsorted > merged_CADDv1.5_offline.vcf && \
        rm -f merged_CADDv1.5_offline_unsorted && \
        python /usr/src/app/CADD_offline_mito_postprocessing.py -c /usr/src/app/whole_mito_SNP_pp2_predictions_sorted.txt -i merged_CADDv1.5_offline.vcf -o merged_CADDv1.5_offline_proper_format.vcf && \
        rm -f merged_CADDv1.5_offline.vcf
    >>>
    output {
        File merged_cadd_output_vcf = "merged_CADDv1.5_offline_proper_format.vcf"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 4
        docker: "quay.io/cmarkello/cadd:latest"
    }
}

task run_cadd_editor {
    input {
        File in_vs_file
        File? in_cadd_output
    }
    command <<<
        set -exu -o pipefail
        
        java -cp /cadd_edit/NewCaddEditor.jar:/cadd_edit/commons-cli-1.4.jar NewCaddEditor --input_vs ~{in_vs_file} --output_cadd ~{in_cadd_output} --output_vs "cadd_editor_output.vs"
    >>>
    output {
        File cadd_editor_output_vs = "cadd_editor_output.vs"
    }
    runtime {
        memory: 100 + " GB"
        cpu: 4
        docker: "quay.io/cmarkello/cadd_editor:latest"
    }
}

