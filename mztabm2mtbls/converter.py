import click
import json
import os
import sys
from typing import List

from metabolights_utils import IsaTableFileReaderResult
from metabolights_utils.isatab import Reader, Writer
from metabolights_utils.models.isa.assay_file import AssayFile
from metabolights_utils.models.isa.assignment_file import AssignmentFile
from metabolights_utils.models.isa.investigation_file import (
    Assay, BaseSection, Factor, Investigation, InvestigationContacts,
    InvestigationPublications, OntologyAnnotation, OntologySourceReference,
    OntologySourceReferences, Person, Protocol, Publication, Study,
    StudyAssays, StudyContacts, StudyFactors, StudyProtocols,
    StudyPublications, ValueTypeAnnotation)
from metabolights_utils.models.isa.samples_file import SamplesFile
from metabolights_utils.models.metabolights.model import MetabolightsStudyModel

from mztabm2mtbls import utils
from mztabm2mtbls.mapper.base_mapper import BaseMapper
from mztabm2mtbls.mapper.metadata.metadata_assay import MetadataAssayMapper
from mztabm2mtbls.mapper.metadata.metadata_base import MetadataBaseMapper
from mztabm2mtbls.mapper.metadata.metadata_contact import MetadataContactMapper
from mztabm2mtbls.mapper.metadata.metadata_cv import MetadataCvMapper
from mztabm2mtbls.mapper.metadata.metadata_database import \
    MetadataDatabaseMapper
from mztabm2mtbls.mapper.metadata.metadata_publication import \
    MetadataPublicationMapper
from mztabm2mtbls.mapper.metadata.metadata_sample import MetadataSampleMapper
from mztabm2mtbls.mapper.metadata.metadata_sample_processing import \
    MetadataSampleProcessingMapper
from mztabm2mtbls.mapper.metadata.metadata_software import \
    MetadataSoftwareMapper
from mztabm2mtbls.mapper.summary.small_molecule_summary import \
    SmallMoleculeSummaryMapper
from mztabm2mtbls.mztab2 import MzTab

mappers: List[BaseMapper] = [
    MetadataBaseMapper(),
    MetadataContactMapper(),
    MetadataPublicationMapper(),
    MetadataCvMapper(),
    MetadataSampleMapper(),
    MetadataSampleProcessingMapper(),
    MetadataSoftwareMapper(),
    MetadataDatabaseMapper(),
    MetadataAssayMapper(),
    SmallMoleculeSummaryMapper()
]
        
@click.command()
@click.option('--input-file', default='test/data/singaporean-plasma-site1.mzTab.json', help='The mzTab-M file in JSON format to convert.')
@click.option('--output_dir', default='output', help='The directory to save the converted files.')
@click.option('--mtbls_accession_number', default='MTBLS1000000', help='The MetaboLights study accession number.')     
def convert(input_file: str, output_dir: str, mtbls_accession_number: str):

    #print disclaimer that we currently do not fully validate neither the mzTab-M file, nor the ISA-Tab files
    print("Please note that the mzTab-M file is not fully validated by this tool. The ISA-Tab files are not validated either at the moment.")

    if input_file.startswith("test/data/singa"):
        print("Using default input file: test/data/singaporean-plasma-site1.mzTab.json")
        print("Please use the --input-file option to specify a custom input file.")
        
    with open(input_file) as f:
        mztab_json_data = json.load(f)
    utils.replace_null_string_with_none(mztab_json_data)
    mztab_model: MzTab = MzTab.model_validate(mztab_json_data)
    utils.modify_mztab_model(mztab_model)
    mtbls_model: MetabolightsStudyModel = utils.create_metabolights_study_model(study_id=mtbls_accession_number)

    for mapper in mappers:
        mapper.update(mztab_model, mtbls_model)
    study_metadata_output_path = os.path.join(output_dir, mtbls_accession_number)
    utils.save_metabolights_study_model(mtbls_model, output_dir=study_metadata_output_path)

if __name__ == "__main__":
    convert()
