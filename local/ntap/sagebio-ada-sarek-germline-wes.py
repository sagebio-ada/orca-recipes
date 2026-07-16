"""https://sagebionetworks.jira.com/browse/EAGER-germline-1
Orca recipe to run the `nf-core/sarek` pipeline for GRCh38 Whole Exome Sequencing (WES) germline variant calling.

Orchestrates four steps:
  1. stage_samplesheet   : Download samplesheet from Synapse and upload to S3
  2. nf-synapse SYNSTAGE : Download FASTQ files from Synapse to S3
  3. nf-core/sarek       : Run germline variant calling
  4. nf-synapse SYNINDEX : Index results back to Synapse

Prerequisites:
  - pip install py-orca
  - AWS credentials configured (profile: tower)
  - SYNAPSE_AUTH_TOKEN set in Tower workspace secrets (not user secret)
    (shared group Synapse token — authenticates the workflow to Synapse;
     manage at https://tower.sagebionetworks.org/orgs/Sage-Bionetworks/workspaces/ntap-add5-project/secrets)

Setup:
  Get your Tower personal access token at tower.sagebionetworks.org → Your tokens,
  then export the following before running:

    export TOWER_ACCESS_TOKEN="<your_tower_token>"
    export TOWER_WORKSPACE="sage-bionetworks/ntap-add5-project"
    export TOWER_API_ENDPOINT="https://tower.sagebionetworks.org/api"
    export NEXTFLOWTOWER_CONNECTION_URI="https://:<your_tower_token>@tower.sagebionetworks.org/api?workspace=sage-bionetworks%2Fntap-add5-project"

Usage:
  # Run all steps (default)
  python sagebio-ada-sarek-germline-wes.py
  python sagebio-ada-sarek-germline-wes.py --run-number 2 # with a specific run number


  # Run individual steps (e.g. stage_samplesheet, synstage, process, synindex)
  python sagebio-ada-sarek-germline-wes.py stage_samplesheet
  python sagebio-ada-sarek-germline-wes.py synstage
  python sagebio-ada-sarek-germline-wes.py process
  python sagebio-ada-sarek-germline-wes.py synindex

  # Run multiple steps (e.g. synstage + process)
  python sagebio-ada-sarek-germline-wes.py synstage process
  python sagebio-ada-sarek-germline-wes.py synstage process --run-number 2 # with a specific run number

  # Rerun with a new version number to preserve previous S3 outputs
  # Default is --run-number 1; increment each time you need a clean rerun
  python sagebio-ada-sarek-germline-wes.py process --run-number 2
  python sagebio-ada-sarek-germline-wes.py synindex --run-number 2
"""
import asyncio
import argparse
from dataclasses import dataclass

import boto3
from orca.services.nextflowtower import NextflowTowerOps
from orca.services.nextflowtower.models import LaunchInfo
from synapseclient import Synapse

session = boto3.Session(profile_name="tower")
s3 = session.client("s3")


async def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('step', nargs='*', default='all', help='Processing step (Default: all)')
    parser.add_argument('--run-number', type=int, default=1, help='Run version number (Default: 1). Increment to preserve previous outputs.')
    args = parser.parse_args()

    ops = NextflowTowerOps()
    datasets = generate_datasets(run_number=args.run_number)
    runs = [run_workflows(ops, dataset, args.step) for dataset in datasets]
    statuses = await asyncio.gather(*runs)
    print(statuses)


@dataclass
class Dataset:
    id: str
    """The synapse id for the samplesheet."""

    samplesheet: str
    """The name of the samplesheet to run."""

    synapse_id_for_output: str
    """The synapse id for the output folder, this is where the output will be uploaded to."""

    bucket_name: str
    """The name of the bucket to stage the samplesheet in."""

    staging_key: str
    """The key in the S3 bucket where this workflow is going to run."""

    institution: str
    """The institution that generated the samples ('JH' or 'WU'). Determines the BED file used."""

    run_number: int = 1
    """Run version number. Passed from CLI --run-number; increment to preserve previous outputs."""

    @property
    def intervals(self) -> str:
        """The S3 uri for the BED file, determined by institution."""
        if self.institution == "JH":
            return BED_JH
        elif self.institution == "WU":
            return BED_WU

    @property
    def samplesheet_location(self) -> str:
        """The location where the unstaged samplesheet is located."""
        return f"{self.samplesheet_location_prefix}{self.samplesheet}"

    @property
    def samplesheet_to_stage_key(self) -> str:
        """The key in the S3 bucket where the samplesheet is going to be staged."""
        return f"{self.staging_key}to_stage/{self.samplesheet}"

    @property
    def staged_samplesheet_location(self) -> str:
        """The S3 uri where the samplesheet is staged."""
        return f"{self.staging_location}synstage_{self.id}/{self.samplesheet}"

    @property
    def staging_location(self) -> str:
        """The S3 uri where the workflow is going to be run."""
        return f"s3://{self.bucket_name}/{self.staging_key}"

    @property
    def samplesheet_location_prefix(self) -> str:
        """The S3 uri where the unstaged samplesheet is located."""
        return f"s3://{self.bucket_name}/{self.staging_key}to_stage/"

    @property
    def output_directory(self) -> str:
        """The S3 uri where the output is going to be uploaded to. The is used as the
        input for the synindex workflow."""
        return f"s3://{self.bucket_name}/outputs/sarek_GRCh38_{self.id}_{self.run_number}/"

    @property
    def synstage_run_name(self) -> str:
        """The name of the synstage run."""
        return f"synstage_{self.id}"

    @property
    def sarek_run_name(self) -> str:
        """The name of the sarek run."""
        return f"sarek_GRCh38_{self.id}_{self.run_number}"

    @property
    def synindex_run_name(self) -> str:
        """The name of the synindex run."""
        return f"synindex_{self.id}_{self.run_number}"


# BED files for exome seq data from JHU NF1 repository - different batches/institutions
BED_JH = "s3://ntap-add5-project-tower-bucket/reference/Baits_BED_Files_AgilentV6_REVISED_S07604514_ALLBED_merged_020816_withChr_GRCh38_sorted.bed"
BED_WU = "s3://ntap-add5-project-tower-bucket/reference/xgen-exome-research-panel-v2-probes-hg3862a5791532796e2eaa53ff00001c1b3c.bed"


def generate_datasets(run_number: int = 1) -> list[Dataset]:
    """Generate list of datasets.

    Source: https://sagebionetworks.jira.com/browse/WORKFLOWS-538

    Samplesheets: https://www.synapse.org/Synapse:syn74378396
    """
    return [
        Dataset(
            # JH_batch1 tumors, single-lane
            # 26 samples, 26 rows, status = 1 (tumor)
            id="syn74378522",
            samplesheet="sarek_samplesheet_JH_batch1_tumor_germline.csv",
            staging_key="samplesheets/Sarek_Process/EAGER-germline/",
            bucket_name="ntap-add5-project-tower-bucket",
            synapse_id_for_output="syn74391336",
            institution="JH",
            run_number=run_number,
        ),
        Dataset(
            # WU_batch1 tumors, single-lane
            # 26 samples, 26 rows, status = 1 (tumor)
            id="syn74722658",
            samplesheet="sarek_samplesheet_WU_batch1_tumor_germline.csv",
            staging_key="samplesheets/Sarek_Process/EAGER-germline/",
            bucket_name="ntap-add5-project-tower-bucket",
            synapse_id_for_output="syn74530464",
            institution="WU",
            run_number=run_number,
        ),
        Dataset(
            # WU_batch2 tumors, single-lane
            # 20 samples, 20 rows, status = 1 (tumor)
            id="syn74378526",
            samplesheet="sarek_samplesheet_WU_batch2_tumor_germline.csv",
            staging_key="samplesheets/Sarek_Process/EAGER-germline/",
            bucket_name="ntap-add5-project-tower-bucket",
            synapse_id_for_output="syn74530478",
            institution="WU",
            run_number=run_number,
        ),
        Dataset(
            # WU_batch3 blood normals, 2 lanes each
            # 16 samples, 32 rows, status = 0 (normal)
            id="syn74378521",
            samplesheet="sarek_samplesheet_WU_batch3_normal_germline.csv",
            staging_key="samplesheets/Sarek_Process/EAGER-germline/",
            bucket_name="ntap-add5-project-tower-bucket",
            synapse_id_for_output="syn74530480",
            institution="WU",
            run_number=run_number,
        ),
        Dataset(
            # WU_batch3 tumors, single-lane
            # 31 samples, 31 rows, status = 1 (tumor)
            id="syn74378531",
            samplesheet="sarek_samplesheet_WU_batch3_tumor_germline.csv",
            staging_key="samplesheets/Sarek_Process/EAGER-germline/",
            bucket_name="ntap-add5-project-tower-bucket",
            synapse_id_for_output="syn74530479",
            institution="WU",
            run_number=run_number,
        ),
        Dataset(
            # WU_batch_mismatched tumors, single-lane
            # 2 samples, 2 rows, status = 1 (tumor)
            id="syn74385711",
            samplesheet="sarek_samplesheet_WU_batch_mismatched_tumor_germline.csv",
            staging_key="samplesheets/Sarek_Process/EAGER-germline/",
            bucket_name="ntap-add5-project-tower-bucket",
            synapse_id_for_output="syn74530491",
            institution="WU",
            run_number=run_number,
        ),
    ]


def stage_samplesheet(syn: Synapse, dataset: Dataset) -> None:
    """Download the samplesheet from synapse and upload it to S3 in the location where synstage
    is going to grab the file.

    Arguments:
        syn: The logged in synapse instance
        dataset: The dataset to stage the samplesheet for
    """
    samplesheet_file = syn.get(dataset.id)
    samplesheet_file_path = samplesheet_file.path

    s3.upload_file(
        samplesheet_file_path, dataset.bucket_name, dataset.samplesheet_to_stage_key
    )


def prepare_synstage_info(dataset: Dataset) -> LaunchInfo:
    """Generate LaunchInfo for nf-synstage.

    Arguments:
        dataset: The dataset to stage the samplesheet for

    Returns:
        The Nextflow Tower workflow launch specification for synstage step
    """
    return LaunchInfo(
        run_name=dataset.synstage_run_name,
        pipeline="Sage-Bionetworks-Workflows/nf-synapse",
        revision="main",
        profiles=["sage"],
        params={
            "input": dataset.samplesheet_location,
            "outdir": dataset.staging_location,
            "entry": "synstage",
        },
        workspace_secrets=["SYNAPSE_AUTH_TOKEN"]  # set as workspace secret (not user secret) in Tower
    )


def prepare_sarek_launch_info(dataset: Dataset) -> LaunchInfo:
    """Generate LaunchInfo for nf-core/sarek workflow run.

    Arguments:
        dataset: The dataset to stage the samplesheet for

    Returns:
        The Nextflow Tower workflow launch specification for sarek processing step
    """
    return LaunchInfo(
        run_name=dataset.sarek_run_name,
        pipeline="nf-core/sarek",
        revision="3.2.2",
        profiles=["sage"],
        params={
            "input": dataset.staged_samplesheet_location,
            "outdir": dataset.output_directory,
            "wes": True,
            "intervals": dataset.intervals,
            "igenomes_base": "s3://sage-igenomes/igenomes",
            "genome": "GATK.GRCh38",
            "tools": "strelka",
        }
    )


def prepare_synindex_launch_info(dataset: Dataset) -> LaunchInfo:
    """Generate LaunchInfo for nf-synindex workflow run.

    Arguments:
        dataset: The dataset to stage the samplesheet for

    Returns:
        The Nextflow Tower workflow launch specification for synindex step
    """
    return LaunchInfo(
        run_name=dataset.synindex_run_name,
        pipeline="Sage-Bionetworks-Workflows/nf-synapse",
        revision="main",
        profiles=["sage"],
        params={
            "s3_prefix": dataset.output_directory,
            "parent_id": dataset.synapse_id_for_output,
            "entry": "synindex",
        },
        workspace_secrets=["SYNAPSE_AUTH_TOKEN"]  # set as workspace secret (not user secret) in Tower
    )


async def run_workflows(ops: NextflowTowerOps, dataset: Dataset, step):
    if 'all' in step or 'stage_samplesheet' in step:
        print('staging samplesheet')
        syn = Synapse()
        syn.login()
        stage_samplesheet(syn, dataset)

    if 'all' in step or 'synstage' in step:
        print('starting synstage')
        synstage_info = prepare_synstage_info(dataset)
        synstage_run_id = ops.launch_workflow(synstage_info, "spot", ignore_previous_runs=True)
        status = await ops.monitor_workflow(run_id=synstage_run_id, wait_time=60 * 2)
        print(status)

    if 'all' in step or 'process' in step:
        print('starting data processing pipeline')
        sarek_info = prepare_sarek_launch_info(dataset)
        sarek_run_id = ops.launch_workflow(sarek_info, "spot", ignore_previous_runs=True)
        status = await ops.monitor_workflow(run_id=sarek_run_id, wait_time=60 * 2)
        print(status)

    if 'all' in step or 'synindex' in step:
        print('starting synindex')
        synindex_info = prepare_synindex_launch_info(dataset)
        synindex_run_id = ops.launch_workflow(synindex_info, "spot", ignore_previous_runs=True)
        status = await ops.monitor_workflow(run_id=synindex_run_id, wait_time=60 * 2)
        print(status)


if __name__ == "__main__":
    asyncio.run(main())
