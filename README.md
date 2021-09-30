UDPAnalysisBMTB
---------------
UDP Candidate Analysis scripts and WDL workflows for detecting causal genetic variation

## NIH Biowulf Usage
#### Setup the main working directory
```
cd /data/$USER
```
Launch an interactive session on Biowulf and load requisite Biowulf modules:
```
sinteractive
module load cromwell/40 git python/3.6
```
Clone the github repo and create a work directory for running the wdl workflow:
```
BMTB_WDL_DIR="/data/$USER/bmtb_wdl_run/wdl_tools"
mkdir -p ${BMTB_WDL_DIR} && cd ${BMTB_WDL_DIR}
git clone https://github.com/cmarkello/bmtb_wdl.git
```
Create workflow inputs directory, download workflow inputs, and setup miniwdl virtual environment:
```
WORKFLOW_INPUT_DIR="/data/$USER/bmtb_wdl_run/workflow_inputs"
mkdir -p ${WORKFLOW_INPUT_DIR}
chmod 2770 ${WORKFLOW_INPUT_DIR}
cd ${BMTB_WDL_DIR}
git clone https://github.com/cmarkello/miniwdl.git
virtualenv miniwdl_venv
source miniwdl_venv/bin/activate
pip3 install ./miniwdl
deactivate
exit
```
#### Running a pedigree sample
Setup workflow directory if it doesn't already exist
```
COHORT_WORKFLOW_DIR="/data/$USER/bmtb_wdl_run/cohort_workdir"
mkdir -p ${COHORT_WORKFLOW_DIR} cd ${COHORT_WORKFLOW_DIR}
```
Launch an interactive node with at least 100GB of memory
```
sinteractive -c6 --mem=100gb
module load cromwell/40 python/3.6
source /data/$USER/bmtb_wdl_run/wdl_tools/miniwdl_venv/bin/activate
cd /data/$USER/bmtb_wdl_run/cohort_workdir
```
Set requisite workflow inputs
```
COHORT_VCF_FILE_PATH="/PATH/TO/JOINT/CALLED/VARSIFTER/FILE.vcf.gz"
MATERNAL_INPUT_BAM_FILE="/PATH/TO/MATERNAL/BAM/FILE.bam"
MATERNAL_INPUT_BAM_FILE_INDEX="/PATH/TO/MATERNAL/BAM/FILE/INDEX.bam.bai"
PATERNAL_INPUT_BAM_FILE="/PATH/TO/PATERNAL/BAM/FILE.bam"
PATERNAL_INPUT_BAM_FILE_INDEX="/PATH/TO/PATERNAL/BAM/FILE/INDEX.bam.bai"
PROBAND_BAM_FILE="/PATH/TO/PROBAND/BAM/FILE.bam"
PROBAND_BAM_FILE_INDEX="/PATH/TO/PROBAND/BAM/FILE/INDEX.bam.bai"
SIBLING_1_BAM_FILE="/PATH/TO/SIBLING/BAM/FILE.bam"
SIBLING_1_BAM_FILE_INDEX="/PATH/TO/SIBLING/BAM/FILE/INDEX.bam.bai"
SAMPLE_NAME_MATERNAL="MATERNAL_UDP_SAMPLE_NAME"
SAMPLE_NAME_PATERNAL="PATERNAL_UDP_SAMPLE_NAME"
SAMPLE_NAME_PROBAND="PROBAND_UDP_SAMPLE_NAME"
SAMPLE_NAME_SIBLING_1="SIBLING_UDP_SAMPLE_NAME"
GENDER_PROBAND=0                                  ## (0 = male, 1 = female)
GENDER_SIBLING_1=1
AFFECTED_PROBAND=0                                ## (1 = unaffected, 0 = affected)
AFFECTED_SIBLING=1
```
Run the wdl workflow
NOTE: Depending on the number of copies of the `SIBLING_BAM_FILE_LIST`, `SIBLING_BAM_FILE_INDEX_LIST`, `SAMPLE_NAME_SIBLING_LIST`, `GENDER_SIBLING_LIST`, and `AFFECTED_SIBLING_LIST` arguments depends on the number of siblings in the pedigree dataset.
```
miniwdl cromwell /data/$USER/bmtb_wdl_run/wdl_tools/bmtb_wdl/wdl/bmtb.wdl \
INPUT_VCF_FILE=${COHORT_VCF_FILE_PATH} \
MATERNAL_INPUT_BAM_FILE=${MATERNAL_INPUT_BAM_FILE} \
MATERNAL_INPUT_BAM_FILE_INDEX=${MATERNAL_INPUT_BAM_FILE_INDEX} \
PATERNAL_INPUT_BAM_FILE=${PATERNAL_INPUT_BAM_FILE} \
PATERNAL_INPUT_BAM_FILE_INDEX=${PATERNAL_INPUT_BAM_FILE_INDEX} \
SIBLING_BAM_FILE_LIST=${PROBAND_BAM_FILE} \
SIBLING_BAM_FILE_LIST=${SIBLING_1_BAM_FILE} \
SIBLING_BAM_FILE_INDEX_LIST=${PROBAND_BAM_FILE_INDEX} \
SIBLING_BAM_FILE_INDEX_LIST=${SIBLING_1_BAM_FILE_INDEX} \
SAMPLE_NAME_MATERNAL=${SAMPLE_NAME_MATERNAL} \
SAMPLE_NAME_PATERNAL=${SAMPLE_NAME_PATERNAL} \
SAMPLE_NAME_SIBLING_LIST=${SAMPLE_NAME_PROBAND} \
SAMPLE_NAME_SIBLING_LIST=${SAMPLE_NAME_SIBLING_1} \
GENDER_SIBLING_LIST=${GENDER_PROBAND} \
GENDER_SIBLING_LIST=${GENDER_SIBLING_1} \
AFFECTED_SIBLING_LIST=${AFFECTED_PROBAND} \
AFFECTED_SIBLING_LIST=${AFFECTED_SIBLING} \
-c /data/$USER/bmtb_wdl_run/wdl_tools/bmtb_wdl/wdl/custom_biowulf_cromwell_singularity.conf \
-d bmtb.final_outputs
```
Output .tar file containing all Black Magic Toolbox analysis data can be found in the following directory
```
/data/$USER/bmtb_wdl_run/cohort_workdir/outputs/bmtb.final_outputs/${SAMPLE_NAME_PROBAND}_BlackBox_Output.tar.gz
```

