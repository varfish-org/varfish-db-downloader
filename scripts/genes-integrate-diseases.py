"""Create integrated gene to disease map."""

import csv
import enum
import gzip
import json
import os
import pickle
import re
import sys
from itertools import chain
from typing import Dict, List, Optional, Set, Tuple

import pronto
from loguru import logger
from pydantic import BaseModel, ConfigDict, Field

#: Path to HGNC:ID to OMIM ID mapping file.
MIM2GENE_MEDGEN_PATH = os.environ["MIM2GENE_MEDGEN_PATH"]
#: Path to genes xlink file.
GENES_XLINK_PATH = os.environ["GENES_XLINK_PATH"]
#: Path to HPOA file.
HPOA_PATH = os.environ["HPOA_PATH"]
#: Path to the ORPHA JSONL file.
ORPHA_JSONL_PATH = os.environ["ORPHA_JSONL_PATH"]
#: Path to the PanelApp JSONL file.
PANELAPP_JSONL_PATH = os.environ["PANELAPP_JSONL_PATH"]
#: Path to the MONDO OBO file.
MONDO_OBO_PATH = os.environ["MONDO_OBO_PATH"]
#: Path to the MONDO unmapped OMIM file.
MONDO_UNMAPPED_OMIM_PATH = os.environ["MONDO_UNMAPPED_OMIM_PATH"]
#: Path to CTD file.
CTD_PATH = os.environ["CTD_PATH"]
#: DO OMIM "ummapped" file.
DO_OMIM_UNMAPPED_PATH = os.environ["DO_OMIM_UNMAPPED_PATH"]
#: DO OMIMinDO file.
DO_OMIM_INDO_PATH = os.environ["DO_OMIM_INDO_PATH"]
#: DO legacy OMIM import file.
DO_OMIM_IMPORT_PATH = os.environ["DO_OMIM_IMPORT_PATH"]

#: Development mode => pickling of data
DEV_MODE = os.environ.get("DEV_MODE", "0") == "1"


class GeneXlink(BaseModel):
    """Cross-reference of gene."""

    #: HGNC ID.
    hgnc_id: str
    #: ENSEMBL gene ID.
    ensembl_gene_id: Optional[str]
    #: Entrez gene ID.
    entrez_id: Optional[str]
    #: Gene symbol
    gene_symbol: str
    #: OMIM IDs
    omim_ids: List[str]


def parse_gene_xlink(path: str) -> List[GeneXlink]:
    """Parse gene cross-reference file."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for key in ("ensembl_gene_id", "entrez_id"):
                if row[key] == "":
                    row[key] = None
            if not row["omim_ids"]:
                row["omim_ids"] = []
            else:
                row["omim_ids"] = row["omim_ids"].split(",")
            rows.append(GeneXlink(**row))
    return rows


class Mim2geneMedgen(BaseModel):
    """Mapping from OMIM ID to gene."""

    model_config = ConfigDict(frozen=True)

    #: OMIM ID.
    omim_id: str
    #: Gene HGNC ID.
    hgnc_id: str

    def __lt__(self, other: "Mim2geneMedgen") -> bool:
        """Compare mim2gene_medgen entries."""
        return (self.hgnc_id, self.omim_id) < (other.hgnc_id, other.omim_id)


def parse_mim2gene_medgen(
    path: str, xlink_by_entrez_id: Dict[str, GeneXlink]
) -> List[Mim2geneMedgen]:
    """Parse mim2gene_medgen file."""
    result = []
    with open(path, "rt") as inputf:
        reader = csv.DictReader(inputf, delimiter="\t")
        for row in reader:
            if row["GeneID"] == "-" or row["type"] != "phenotype":
                continue
            if row["GeneID"] not in xlink_by_entrez_id:
                logger.warning(f"Skipping NCBI Gene ID = {row['GeneID']}")
                continue
            result.append(
                Mim2geneMedgen(
                    omim_id=f"OMIM:{row['#MIM number']}",
                    hgnc_id=xlink_by_entrez_id[row["GeneID"]].hgnc_id,
                )
            )
    return list(sorted(set(result)))


@enum.unique
class Evidence(enum.Enum):
    """Evidence from HPOA record"""

    #: Has been inferred by parsing OMIM.
    INFERRED_FROM_ELECTRONIC_ANNOTATION = "IEA"
    #: Extracted from articles in medical literature.
    PUBLISHED_CLINICAL_STUDY = "PCS"
    #: Extracted from knowledge bases such as OMIM and Orphanet.
    TRACEABLE_AUTHOR_STATEMENT = "TAS"


@enum.unique
class Sex(enum.Enum):
    """Enumeration for sex."""

    #: Male
    MALE = "MALE"
    #: Female
    FEMALE = "FEMALE"


@enum.unique
class Aspect(enum.Enum):
    """Enumeration for aspect."""

    #: Phenotype
    PHENOTYPE = "P"
    #: Mode of inheritance
    MODE_OF_INHERITANCE = "I"
    #: Clinical course
    CLINICAL_COURSE = "C"
    #: Modifier
    CLINICAL_MODIFIER = "M"


class HpoaEntry(BaseModel):
    """Disease-phenotype association from HPOA file."""

    model_config = ConfigDict(use_enum_values=True)

    #: Term database ID.
    database_id: str
    #: Term database name.
    disease_name: str
    #: Whether is a positive or negative association.
    qualifier: bool
    #: HPO ID.
    hpo_id: str
    #: Evidence for this association.
    reference: List[str]
    #: Evidence for this association.
    evidence: Evidence
    #: HPO term for onset.
    onset: str
    #: Frequency annotation.
    frequency: str
    #: Sex link, if any.
    sex: Optional[Sex]
    #: HPO terms of modifiers.
    modifier: List[str]
    #: Aspect of this association.
    aspect: Aspect
    #: Biocuration notes.
    biocuration: List[str]


def parse_hpoa(path: str) -> List[HpoaEntry]:
    """Parse HPOA file."""
    rows = []
    with open(path) as f:
        for i in range(4):
            next(f)
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row["qualifier"] = row["qualifier"] != "NOT"
            row["reference"] = row["reference"].split(";")
            row["evidence"] = Evidence(row["evidence"])
            if row["sex"]:
                row["sex"] = Sex(row["sex"].upper())
            else:
                row["sex"] = None
            row["modifier"] = row["modifier"].split(";")
            row["aspect"] = Aspect(row["aspect"])
            row["biocuration"] = row["biocuration"].split(";")
            rows.append(HpoaEntry(**row))
    return rows


class DiseaseNameLabel(BaseModel):
    """Association of a diasease database ID to its name."""

    #: Disease database ID.
    database_id: str
    #: Name in the database.
    disease_name: str


def hpoa_entries_to_disease_label_map(hpoa_entries: List[HpoaEntry]) -> List[DiseaseNameLabel]:
    """Create a map of disease database ID to disease name."""
    disease_label_map = {}
    for entry in hpoa_entries:
        disease_label_map[entry.database_id] = DiseaseNameLabel(
            database_id=entry.database_id, disease_name=entry.disease_name
        )
    return disease_label_map


@enum.unique
class DisorderMappingRelation(enum.Enum):
    """Characterisation of disorder mapping relation."""

    #: Exact mapping
    EXACT = "E"
    #: ORPHA code narrower term maps to a broader term
    NTBT = "NTBT"
    #: ORPHA code broader term maps to a narrower term
    BTNT = "BTNT"
    #: ORPHA code narrower term maps to a broader term because of an exact
    #: match with an synonym in the target terminology.
    NTBT_E = "NTBT/E"
    #: ORPHA code broader term maps to a narrower term because of an exact
    #: match with an synonym in the target terminology.
    BTNT_E = "BTNT/E"
    #: Not yet decied / unable to decide.
    ND = "ND"


class OrphaDisorderMapping(BaseModel):
    """Mapping into a external/target disorder terminology."""

    model_config = ConfigDict(use_enum_values=True)

    #: Identifier of source terminology with prefix of the terminology.
    database_id: str
    #: The relationship between ORPHA and the target terminology.
    relation: DisorderMappingRelation


class OrphaDisorder(BaseModel):
    """Description of the ORPHA disorder."""

    #: The "ORPHA:<ID>" ID.
    database_id: str
    #: The name of the disorder.
    database_name: str
    #: Definition, if provided.
    definition: Optional[str]
    #: Synonym list.
    synonyms: List[str]
    #: Mappings in external terminologies.
    mappings: List[OrphaDisorderMapping]


@enum.unique
class OrphaStatus(enum.Enum):
    """ORPHA assessment status."""

    #: Assessed
    ASSESSED = "Assessed"
    #: Not yet assessed
    NOT_YET_ASSESSED = "Not yet assessed"


@enum.unique
class OrphaAssocationType(enum.Enum):
    """ORPHA association type."""

    ROLE = "Role in the phenotype of"
    GERMLINE_LOF = "Disease-causing germline mutation(s) (loss of function) in"
    GERMLINE_GOF = "Disease-causing germline mutation(s) (gain of function) in"
    GERMLINE = "Disease-causing germline mutation(s) in"
    SOMATIC = "Disease-causing somatic mutation(s) in"
    MODIFYING_GERMLINE = "Modifying germline mutation in"
    FUSION_GENE = "Part of a fusion gene in"
    SUSCEPTIBILITY = "Major susceptibility factor in"
    CANDIDATE_GENE = "Candidate gene tested in"
    BIOMARKER = "Biomarker tested in"


class OrphaGeneMapping(BaseModel):
    """Mapping from an ORPHA disorder to a gene."""

    model_config = ConfigDict(use_enum_values=True)

    #: ORPHA ID.
    orpha_id: str
    #: Gene HGNC id.
    hgnc_id: str
    #: Assessment status.
    status: OrphaStatus
    #: The type of the assessment.
    association_type: OrphaAssocationType


def parse_orpha_jsonl(
    path: str, xlink_by_gene_symbol: Dict[str, GeneXlink]
) -> Tuple[List[OrphaDisorder], List[OrphaGeneMapping]]:
    """Parse ORPHA JSONL file."""
    disorders = []
    gene_mappings = []
    with open(path) as inputf:
        for line in inputf:
            record = json.loads(line)
            # Extract disease details and cross-references to other terminologies.
            disorder_mappings = []
            for mapping in record["cross_references"].get("ExternalReference") or []:
                disorder_mappings.append(
                    OrphaDisorderMapping(
                        database_id=f"{mapping['Source']}:{mapping['Reference']}",
                        relation=DisorderMappingRelation(
                            mapping["DisorderMappingRelation"].split(" ")[0]
                        ),
                    )
                )
            summary_information = record["cross_references"].get("SummaryInformation") or []
            disorder = OrphaDisorder(
                database_id=record["orpha_id"],
                database_name=record["cross_references"]["Preferred term"],
                definition=summary_information[0].get("Definition")
                if summary_information
                else None,
                synonyms=record["cross_references"].get("Synonym") or [],
                mappings=disorder_mappings,
            )
            disorders.append(disorder)
            # Extract gene mappings.
            for assoc in (record.get("disease_genes") or {}).get("DisorderGeneAssociation") or []:
                hgnc_id = None
                for ext_ref in assoc["Gene"]["ExternalReference"] or []:
                    if ext_ref["Source"] == "HGNC":
                        hgnc_id = f"HGNC:{ext_ref['Reference']}"
                if hgnc_id is None:
                    if assoc["Gene"]["Symbol"] not in xlink_by_gene_symbol:
                        logger.info(f"Skipping GeneType = {assoc['Gene']['GeneType']}")
                        continue
                    else:
                        hgnc_id = xlink_by_gene_symbol[assoc["Gene"]["Symbol"]].hgnc_id
                gene_mappings.append(
                    OrphaGeneMapping(
                        orpha_id=record["orpha_id"],
                        hgnc_id=hgnc_id,
                        status=OrphaStatus(assoc["DisorderGeneAssociationStatus"]),
                        association_type=OrphaAssocationType(assoc["DisorderGeneAssociationType"]),
                    )
                )
    return disorders, gene_mappings


@enum.unique
class PanelappConfidence(enum.Enum):
    """Confidence level for PanelApp."""

    #: Green
    GREEN = "GREEN"
    #: Amber
    AMBER = "AMBER"
    #: Red
    RED = "RED"
    #: None
    NONE = "NONE"

    def __lt__(self, other: "PanelappConfidence") -> bool:
        """Compare confidence levels.

        GREEN is smallest for sorting ergnonomics.
        """
        mapping = {
            PanelappConfidence.GREEN: 0,
            PanelappConfidence.AMBER: 1,
            PanelappConfidence.RED: 2,
            PanelappConfidence.NONE: 3,
        }
        return mapping[self] < mapping[other]


@enum.unique
class PanelappEntityType(enum.Enum):
    """Entity type for PanelApp."""

    #: Gene
    GENE = "GENE"
    #: Region
    REGION = "REGION"
    #: STR
    STR = "STR"


class PanelappPanel(BaseModel):
    """Information about a PanelApp panel."""

    #: Numeric ID
    id: int
    #: Name of the panel.
    name: str
    #: Version of the panel.
    version: str


class PanelappAssociation(BaseModel):
    """A gene-to-disease association from PanelApp."""

    model_config = ConfigDict(use_enum_values=True)

    #: Gene HGNC ID.
    hgnc_id: str
    #: Confidence level.
    confidence_level: PanelappConfidence
    #: Entity type.
    entity_type: PanelappEntityType
    #: Mode of inheritance.
    mode_of_inheritance: Optional[str]
    #: The phenotypes as list of strings.
    phenotypes: List[str]
    #: Information about the panel that is the source.
    panel: PanelappPanel


def parse_panelapp_jsonl(
    path: str, xlink_by_gene_symbol: Dict[str, GeneXlink]
) -> List[PanelappAssociation]:
    entity_type_mapping = {
        "gene": PanelappEntityType.GENE.value,
        "region": PanelappEntityType.REGION.value,
        "str": PanelappEntityType.STR.value,
    }
    confidence_level_mapping = {
        "3": PanelappConfidence.GREEN.value,
        "2": PanelappConfidence.AMBER.value,
        "1": PanelappConfidence.RED.value,
        "0": PanelappConfidence.NONE.value,
    }

    result = []
    with open(path, "rt") as inputf:
        for line in inputf:
            record = json.loads(line)
            entity_type = PanelappEntityType(entity_type_mapping[record["entity_type"]])
            if entity_type not in [PanelappEntityType.GENE, PanelappEntityType.STR]:
                continue
            hgnc_id = record["gene_data"].get("hgnc_id")
            if not hgnc_id:
                hgnc_id = getattr(
                    xlink_by_gene_symbol.get(record["gene_data"]["gene_symbol"]), "hgnc_id", None
                )
            if not hgnc_id:
                logger.warn(
                    f"Skipping record without HGNC ID and unmappable gene symbol: {record['gene_data']['gene_symbol']}"
                )
                continue
            record["confidence_level"] = confidence_level_mapping[record["confidence_level"]]
            assoc = PanelappAssociation(
                hgnc_id=hgnc_id,
                confidence_level=PanelappConfidence(record["confidence_level"]),
                mode_of_inheritance=record["mode_of_inheritance"],
                entity_type=entity_type,
                phenotypes=record["phenotypes"],
                panel=PanelappPanel(
                    id=record["panel"]["id"],
                    name=record["panel"]["name"],
                    version=record["panel"]["version"],
                ),
            )
            if assoc.confidence_level != PanelappConfidence.NONE:
                result.append(assoc)
    return result


@enum.unique
class MondoDiseaseRelation(enum.Enum):
    """Characterisation of disorder mapping relation."""

    #: Exact.
    EXACT = "EXACT"
    #: Related.
    RELATED = "RELATED"
    #: Broad.
    BROAD = "BROAD"
    #: Narrow.
    NARROW = "NARROW"
    #: Abbreviation.
    ABBREVIATION = "ABBREVIATION"
    #: Ambiguous.
    AMBIGUOUS = "AMBIGUOUS"
    #: Deprecated.
    DEPRECATED = "DEPRECATED"
    #: Dubious.
    DUBIOUS = "DUBIOUS"
    #: Exclude.
    EXCLUDE = "EXCLUDE"
    #: Misspelling.
    MISSPELLING = "MISSPELLING"
    #: Non-human.
    NON_HUMAN = "NON-HUMAN"
    #: UK spelling synonym.
    OMO_0003005 = "OMO:0003005"


class MondoSynonym(BaseModel):
    """Mapping into another terminology."""

    model_config = ConfigDict(use_enum_values=True)

    #: Identifier of source terminology with prefix of the terminology.
    database_id: str
    #: The relationship between ORPHA and the target terminology.
    relation: List[MondoDiseaseRelation]


class MondoDisease(BaseModel):
    """Relevant information from a MONDO disease."""

    #: The MONDO term ID.
    mondo_id: str
    #: The name of the disease.
    name: str
    #: A list of synonyms.
    synonyms: List[MondoSynonym]


def parse_mondo_obo(path: str) -> List[MondoDisease]:
    """Parse MONDO OBO file."""
    result = []
    mondo = pronto.Ontology(path)
    for term in mondo.terms():
        synonyms = []
        for synonym in term.synonyms:
            for xref in synonym.xrefs:
                synonyms.append(
                    MondoSynonym(
                        database_id=xref.id,
                        relation=list(map(MondoDiseaseRelation, synonym.scope.split(" "))),
                    )
                )
        result.append(MondoDisease(mondo_id=term.id, name=term.name, synonyms=synonyms))
    return result


class MondoUnmappedOmim(BaseModel):
    """Model for the MONDO unmapped OMIM file."""

    subject_id: str
    subject_label: str


def parse_mondo_unmapped_omim(path: str) -> List[MondoUnmappedOmim]:
    """Parse MONDO unmapped OMIM file."""
    result = []
    with open(path) as inputf:
        for line in inputf:
            if line.startswith("subject_id"):
                pass
            record = line.rstrip("\n").split("\t")
            result.append(MondoUnmappedOmim(subject_id=record[0], subject_label=record[1]))
    return result


@enum.unique
class ConfidenceLevel(enum.Enum):
    """Confidence level for gene-disease association."""

    #: High confidence.
    #:
    #: This corresponds to a match in OMIM, an "Assessed" entry in Orphanet, or a PanelApp entry
    #: with GREEN color.
    HIGH = "HIGH"
    #: Medium confidence.
    #:
    #: This corresponds to a PanelApp entry with AMBER color.
    MEDIUM = "MEDIUM"
    #: Low confidence.
    #:
    #: This corresponds to a "Not yet assessed" entry in Orphanet or a PanelApp entry with RED
    #: color.
    LOW = "LOW"

    def __lt__(self, other: "ConfidenceLevel") -> bool:
        """Compare confidence levels.

        HIGH is smallest for sorting ergnonomics.
        """
        mapping = {
            ConfidenceLevel.HIGH: 3,
            ConfidenceLevel.MEDIUM: 2,
            ConfidenceLevel.LOW: 1,
        }
        return mapping[self] < mapping[other]


@enum.unique
class GeneDiseaseAssociationSource(enum.Enum):
    """Enumeration for gene-disease association source."""

    #: OMIM.
    OMIM = "OMIM"
    #: Orphanet.
    ORPHANET = "ORPHANET"
    #: PanelApp.
    PANELAPP = "PANELAPP"

    def __lt__(self, other: "ConfidenceLevel") -> bool:
        """Compare confidence levels."""
        return self.value < other.value


class GeneDiseaseAssociationEntry(BaseModel):
    """A source gene-disease association entry."""

    #: The source.
    source: GeneDiseaseAssociationSource
    #: The confidence level.
    confidence: ConfidenceLevel


class LabeledDisorder(BaseModel):
    """Disorder identifier with a label."""

    model_config = ConfigDict(frozen=True)

    #: ID of the disorder.
    term_id: str
    #: Label of the disorder.
    title: Optional[str]

    def __lt__(self, other: "LabeledDisorder") -> bool:
        """Compare labeled disorders."""
        return self.term_id < other.term_id


class GeneDiseaseAssociation(BaseModel):
    """Association of a gene to a disease."""

    model_config = ConfigDict(frozen=True)

    #: Gene HGNC ID.
    hgnc_id: str
    #: List of primary disease identifiers that diseases came from.
    labeled_disorders: Tuple[LabeledDisorder, ...]
    #: List of disease identifiers, used for clustering.
    disease_ids: Tuple[str, ...] = Field(..., exclude=True)
    #: Disease name for display, if any.
    disease_name: Optional[str]
    #: Disease definition from Orphanet, if available.
    disease_definition: Optional[str]
    #: The overall supporting sources.
    sources: Tuple[GeneDiseaseAssociationSource, ...]
    #: The overall confidence level.
    confidence: ConfidenceLevel

    def merge(self, other: "GeneDiseaseAssociation") -> "GeneDiseaseAssociation":
        if self.hgnc_id != other.hgnc_id:
            raise RuntimeError("Cannot merge different genes.")
        labeled_disease_ids = tuple(
            sorted(set(chain(self.labeled_disorders, other.labeled_disorders)))
        )
        disease_ids = tuple(sorted(set(chain(self.disease_ids, other.disease_ids))))
        disease_name = self.disease_name or other.disease_name
        disease_definition = self.disease_definition or other.disease_definition
        sources = tuple(sorted(set(chain(self.sources, other.sources))))
        if self.confidence < other.confidence:
            confidence = other.confidence
        else:
            confidence = self.confidence
        return GeneDiseaseAssociation(
            hgnc_id=self.hgnc_id,
            labeled_disorders=labeled_disease_ids,
            disease_ids=disease_ids,
            disease_name=disease_name,
            disease_definition=disease_definition,
            sources=sources,
            confidence=confidence,
        )


class CtdDiseaseEntry(BaseModel):
    """Entry in CTD disease database."""

    model_config = ConfigDict(frozen=True)

    disease_name: str
    disease_id: str
    alt_disease_ids: List[str]
    definition: str
    parent_ids: List[str]
    tree_numbers: List[str]
    parent_tree_numbers: List[str]
    synonyms: List[str]
    slim_mappings: List[str]


def parse_ctd_disease_tsv(path: str) -> List[CtdDiseaseEntry]:
    """Parse CTD disease database TSV file."""
    result = []
    with gzip.open(path, "rt") as inputf:
        for line in inputf:
            if line.startswith("#"):
                continue
            record = line.rstrip("\n").split("\t")
            result.append(
                CtdDiseaseEntry(
                    disease_name=record[0],
                    disease_id=record[1],
                    alt_disease_ids=record[2].split("|"),
                    definition=record[3],
                    parent_ids=record[4].split("|"),
                    tree_numbers=record[5].split("|"),
                    parent_tree_numbers=record[6].split("|"),
                    synonyms=record[7].split("|"),
                    slim_mappings=record[8].split("|"),
                )
            )
    return result


class DoLabel(BaseModel):
    """Label for a DO term."""

    #: The OMIM id.
    omim_id: str
    #: The label.
    label: str


def parse_do_omim_unmapped(path: str) -> List[DoLabel]:
    """Parse DO OMIM unmapped file."""
    result = []
    first = True
    with open(path) as inputf:
        for line in inputf:
            record = line.rstrip("\n").split(",")
            if first:
                first = False
            else:
                result.append(DoLabel(omim_id=record[0], label=record[1]))
    return result


class OmimInDo(BaseModel):
    """Entry in OMIMinDO file."""

    model_config = ConfigDict(frozen=True)

    id: str
    label: str
    xrefs: List[str]


def parse_do_omim_in_do(path: str) -> List[OmimInDo]:
    result = []
    with open(path, "rt") as inputf:
        reader = csv.DictReader(inputf, delimiter="\t", escapechar="\\")
        for row in reader:
            row["xrefs"] = [xref for xref in row["xrefs"].split(", ") if xref != "OMIM:genemap2"]
            result.append(OmimInDo(**row))
    return result


def parse_do_omim_import(path: str) -> List[DoLabel]:
    result = []
    omim_import = pronto.Ontology(path)
    for term in omim_import.terms():
        if term.name:
            result.append(DoLabel(omim_id=term.id, label=term.name))
    return result


class GeneConditionAssociation(BaseModel):
    """Diseases and PanelApp information linked to a gene."""

    #: Gene HGNC ID.
    hgnc_id: str
    #: The gene disease association
    disease_associations: Tuple[GeneDiseaseAssociation, ...]
    #: PanelApp panels with related genes.
    panelapp_associations: Tuple[PanelappAssociation, ...]


class ResultContainer(BaseModel):
    """Container for results."""

    #: The results.
    results: Tuple[GeneConditionAssociation, ...]

    def __lt__(self, other: "ResultContainer") -> bool:
        """Compare two result containers."""
        return self.results < other.results


#: Some manual OMIM labels that are on MONDO exclusion lists etc.
MANUAL_OMIM_LABELS = {
    "OMIM:300932": "Thyroxine-binding globulin QTL",
    "OMIM:615018": "Sd(a) polyagglutination syndrome",
    "OMIM:617956": "Beta-glycopyranoside tasting",
    "OMIM:617966": "Low density lipoprotein cholesterol level QTL 7",
    "OMIM:617995": "IMPDH2 enzyme activity, variation in",
    "OMIM:618079": "Low density lipoprotein cholesterol level QTL 8",
    "OMIM:618807": "LPA deficiency, congenital",
    "OMIM:619812": "Blood group, EMM system",
    "OMIM:620116": "Fatty liver disease, protection from",
    "OMIM:620207": "ER blood group system",
}


#: MondoDiseaseRelations that are considered to be "good enough" mappings.
GOOD_MONDO_RELATIONS = (
    MondoDiseaseRelation.EXACT.value,
    MondoDiseaseRelation.NARROW.value,
    MondoDiseaseRelation.BROAD.value,
    MondoDiseaseRelation.RELATED.value,
)
#: DisorderMappingRelation that are considered to be "good enough" mappings.
GOOD_ORPHA_RELATIONS = (
    DisorderMappingRelation.EXACT.value,
    DisorderMappingRelation.BTNT.value,  # ?
    DisorderMappingRelation.BTNT_E.value,
    DisorderMappingRelation.NTBT.value,
    DisorderMappingRelation.NTBT_E.value,  # ?
)


class GeneDiseaseKey(BaseModel):
    """Key for a gene-disease association."""

    model_config = ConfigDict(frozen=True)

    #: Gene HGNC ID.
    hgnc_id: str
    #: Disease database ID.
    disease_id: str


class Integrator:
    """Implementation of the integration algorithm."""

    # data loaded / pickled

    gene_xlinks: List[GeneXlink]
    xlink_by_hgnc_id: Dict[str, GeneXlink]
    xlink_by_gene_symbol: Dict[str, GeneXlink]
    xlink_by_entrez_id: Dict[str, GeneXlink]
    xlink_by_omim_id: Dict[str, GeneXlink]
    mim2gene_medgen: List[Mim2geneMedgen]
    hpoa_entries: List[HpoaEntry]
    disease_label_map: Dict[str, DiseaseNameLabel]
    orpha_disorders: List[OrphaDisorder]
    orpha_mappings: List[OrphaGeneMapping]
    panelapp_associations: List[PanelappAssociation]
    mondo_diseases: List[MondoDisease]
    mondo_unmapped_omim: List[MondoUnmappedOmim]
    ctd_diseases: List[CtdDiseaseEntry]
    do_unmapped_labels: List[DoLabel]
    omim_in_dos: List[OmimInDo]

    # other data

    #: OMIM titles
    omim_titles: Dict[str, str]
    #: Mappings from MONDO to ORPHA
    mondo_to_orpha: Dict[str, Set[str]]
    #: Mappings from MONDO to OMIM
    mondo_to_omim: Dict[str, Set[str]]
    #: Mappings from ORPHA to OMIM
    orpha_to_omim: Dict[str, Set[str]]
    #: Mappings from OMIM to ORPHA
    omim_to_orpha: Dict[str, Set[str]]

    def __init__(self):
        """Initialise the integrator."""
        #: Mapping from `(gene_hgnc_id, disease_id)` to `GeneDiseaseAssociationEntry`.
        self.disease_assocs: Dict[GeneDiseaseKey, GeneDiseaseAssociation] = {}
        #: Mapping from `hgnc_id` to list of `PanelappAssociation`s.
        self.panelapp_assocs: Dict[str, List[PanelappAssociation]] = {}

    def register_disease_assoc(self, assoc: GeneDiseaseAssociation):
        """Register a gene-disease association."""
        found_list = set()
        for disease_id in assoc.disease_ids:
            key = GeneDiseaseKey(hgnc_id=assoc.hgnc_id, disease_id=disease_id)
            if key in self.disease_assocs:
                found_list.add(self.disease_assocs[key])
        if not found_list:
            for disease_id in assoc.disease_ids:
                key = GeneDiseaseKey(hgnc_id=assoc.hgnc_id, disease_id=disease_id)
                self.disease_assocs[key] = assoc
        else:
            if len(found_list) != 1:
                logger.warning(f"Found multiple associations for {assoc.hgnc_id}")
            for found in found_list:
                found = found.merge(assoc)
                for disease_id in assoc.disease_ids:
                    key = GeneDiseaseKey(hgnc_id=assoc.hgnc_id, disease_id=disease_id)
                    self.disease_assocs[key] = found

    def run(self, pickle_path: Optional[str] = None):
        logger.info("Building gene-disease map...")
        self.load_data(pickle_path)
        self.handle_orpha()
        self.handle_mim2gene_medgen()
        self.handle_panelapp()

        conditions_by_hgnc = {
            hgnc_id: GeneConditionAssociation(
                hgnc_id=hgnc_id, disease_associations=[], panelapp_associations=[]
            )
            for hgnc_id in sorted(
                set(
                    chain(
                        (k.hgnc_id for k in self.disease_assocs.keys()), self.panelapp_assocs.keys()
                    )
                )
            )
        }
        for hgnc_id, assocs in self.panelapp_assocs.items():
            conditions_by_hgnc[hgnc_id] = conditions_by_hgnc[hgnc_id].model_copy(
                update={
                    "panelapp_associations": tuple(sorted(assocs, key=lambda a: a.confidence_level))
                }
            )
        seen = set()
        for key, assoc in self.disease_assocs.items():
            disease_assoc = conditions_by_hgnc[key.hgnc_id].disease_associations
            marker = (assoc.hgnc_id, assoc.labeled_disorders)
            if marker not in seen:
                seen.add(marker)
                conditions_by_hgnc[key.hgnc_id] = conditions_by_hgnc[key.hgnc_id].model_copy(
                    update={
                        "disease_associations": tuple(
                            sorted(chain(disease_assoc, [assoc]), key=lambda a: a.confidence)
                        )
                    }
                )
        result = ResultContainer(results=tuple(conditions_by_hgnc.values()))
        for assoc in result.results:
            json.dump(obj=assoc.model_dump(mode="json"), fp=sys.stdout)
            sys.stdout.write("\n")
        logger.info("All done. Have a nice day.")

    def load_data(self, pickle_path: Optional[str] = None):  # noqa: C901
        if pickle_path and os.path.exists(pickle_path):
            logger.info("unpickling...")
            with gzip.open(pickle_path, "rb") as inputf:
                (
                    self.gene_xlinks,
                    self.mim2gene_medgen,
                    self.hpoa_entries,
                    self.disease_label_map,
                    self.orpha_disorders,
                    self.orpha_mappings,
                    self.panelapp_associations,
                    self.mondo_diseases,
                    self.mondo_unmapped_omim,
                    self.ctd_diseases,
                    self.do_unmapped_labels,
                    self.omim_in_dos,
                    self.do_omim_import,
                ) = pickle.load(inputf)
            self.xlink_by_entrez_id = {
                xlink.entrez_id: xlink for xlink in self.gene_xlinks if xlink.entrez_id
            }
            self.xlink_by_gene_symbol = {xlink.gene_symbol: xlink for xlink in self.gene_xlinks}
        else:
            logger.info("Parsing genes xlink file...")
            self.gene_xlinks = parse_gene_xlink(GENES_XLINK_PATH)
            self.xlink_by_entrez_id = {
                xlink.entrez_id: xlink for xlink in self.gene_xlinks if xlink.entrez_id
            }
            logger.info("Parsing mim2gene_medgen file...")
            self.mim2gene_medgen = parse_mim2gene_medgen(
                MIM2GENE_MEDGEN_PATH, self.xlink_by_entrez_id
            )
            logger.info("Parsing HPOA file...")
            self.hpoa_entries = parse_hpoa(HPOA_PATH)
            logger.info("Converting to disease label map...")
            self.disease_label_map = hpoa_entries_to_disease_label_map(self.hpoa_entries)
            # print(list(disease_label_map.items())[:10])
            logger.info("Parsing ORPHA JSONL file...")
            self.xlink_by_gene_symbol = {xlink.gene_symbol: xlink for xlink in self.gene_xlinks}
            self.orpha_disorders, self.orpha_mappings = parse_orpha_jsonl(
                ORPHA_JSONL_PATH, self.xlink_by_gene_symbol
            )
            # print(json.dumps([d.model_dump() for d in disorders[:10]], indent=2))
            logger.info("Parsing PanelApp JSONL file...")
            self.panelapp_associations = parse_panelapp_jsonl(
                PANELAPP_JSONL_PATH, self.xlink_by_gene_symbol
            )

            logger.info("Loading MONDO files...")
            self.mondo_diseases = parse_mondo_obo(MONDO_OBO_PATH)
            self.mondo_unmapped_omim = parse_mondo_unmapped_omim(MONDO_UNMAPPED_OMIM_PATH)
            logger.info("Loading CTD disease file...")
            self.ctd_diseases = parse_ctd_disease_tsv(CTD_PATH)
            logger.info("Loading DO files...")
            self.do_unmapped_labels = parse_do_omim_unmapped(DO_OMIM_UNMAPPED_PATH)
            self.omim_in_dos = parse_do_omim_in_do(DO_OMIM_INDO_PATH)
            self.do_omim_import = parse_do_omim_import(DO_OMIM_IMPORT_PATH)

            if pickle_path:
                with gzip.open(pickle_path, "wb") as outputf:
                    logger.info("Pickling...")
                    pickle.dump(
                        (
                            self.gene_xlinks,
                            self.mim2gene_medgen,
                            self.hpoa_entries,
                            self.disease_label_map,
                            self.orpha_disorders,
                            self.orpha_mappings,
                            self.panelapp_associations,
                            self.mondo_diseases,
                            self.mondo_unmapped_omim,
                            self.ctd_diseases,
                            self.do_unmapped_labels,
                            self.omim_in_dos,
                            self.do_omim_import,
                        ),
                        outputf,
                    )

        #: Build gene mappings.
        self.xlink_by_hgnc_id = {xlink.hgnc_id: xlink for xlink in self.gene_xlinks}
        self.xlink_by_omim_id = {
            f"OMIM:{omim_id}": xlink for xlink in self.gene_xlinks for omim_id in xlink.omim_ids
        }

        # Build mapping from OMIM to title.
        self.omim_titles = {}
        for entry in self.hpoa_entries:
            if entry.database_id.startswith("OMIM:"):
                self.omim_titles[entry.database_id] = entry.disease_name
        for label in self.do_unmapped_labels:
            self.omim_titles.setdefault(label.omim_id, label.label)
        for omim_in_do in self.omim_in_dos:
            for xref in omim_in_do.xrefs:
                if xref.startswith("OMIM:") and xref not in self.omim_titles:
                    self.omim_titles[xref] = omim_in_do.label
        for label in self.do_omim_import:
            if label.omim_id.startswith("OMIM:") and label.omim_id not in self.omim_titles:
                self.omim_titles[label.omim_id] = label.label
        for ctd_disease in self.ctd_diseases:
            for disease_id in chain([ctd_disease.disease_id], ctd_disease.alt_disease_ids):
                if disease_id.startswith("OMIM:") and disease_id not in self.omim_titles:
                    self.omim_titles[disease_id] = ctd_disease.disease_name
        for unmapped_record in self.mondo_unmapped_omim:
            self.omim_titles.setdefault(unmapped_record.subject_id, unmapped_record.subject_label)
        for key, value in MANUAL_OMIM_LABELS.items():
            self.omim_titles.setdefault(key, value)

        # Build ID mappings from MONDO to OMIM and ORPHA.
        self.mondo_to_omim = {}
        self.mondo_to_orpha = {}
        for relation in GOOD_MONDO_RELATIONS:
            for mondo_disease in self.mondo_diseases:
                for synonym in mondo_disease.synonyms:
                    if (
                        synonym.database_id.startswith("OMIM:")
                        and synonym.database_id != "OMIM:genemap2"
                        and relation in synonym.relation
                        and synonym.database_id not in self.mondo_to_omim
                    ):
                        self.mondo_to_omim.setdefault(mondo_disease.mondo_id, set()).add(
                            synonym.database_id
                        )
                    if (
                        synonym.database_id.startswith("Orphanet:")
                        and relation in synonym.relation
                        and synonym.database_id not in self.mondo_to_orpha
                    ):
                        self.mondo_to_orpha.setdefault(mondo_disease.mondo_id, set()).add(
                            synonym.database_id.replace("Orphanet:", "ORPHA:").replace(
                                "-definition", ""
                            )
                        )
        # Build ID mappings between ORPHA and OMIM.
        self.orpha_to_omim = {}
        self.omim_to_orpha = {}
        for orpha_disorder in self.orpha_disorders:
            orpha_id = orpha_disorder.database_id
            for mapping in orpha_disorder.mappings:
                if (
                    mapping.database_id.startswith("OMIM:")
                    and mapping.relation in GOOD_ORPHA_RELATIONS
                ):
                    omim_id = mapping.database_id
                    self.orpha_to_omim.setdefault(orpha_id, set()).add(omim_id)
                    self.omim_to_orpha.setdefault(omim_id, set()).add(orpha_id)

    def handle_mim2gene_medgen(self):
        """Process the mim2gene_medgen data."""
        logger.info("Processing mim2gene_medgen...")
        # Our primary source for OMIM names is the HPOA table.
        hpoa_by_disease_id = {
            entry.database_id: entry for entry in self.hpoa_entries if entry.database_id
        }
        # Our second source for OMIM names are files from DO.
        mondo_by_omim_id: Dict[str, List[MondoDisease]] = {}
        for relation in GOOD_MONDO_RELATIONS:
            for mondo_disease in self.mondo_diseases:
                for synonym in mondo_disease.synonyms:
                    if (
                        synonym.database_id.startswith("OMIM:")
                        and relation in synonym.relation
                        and synonym.database_id not in mondo_by_omim_id
                    ):
                        mondo_by_omim_id.setdefault(synonym.database_id, []).append(mondo_disease)
        # do_by_disease_id: Dict[str, str] = {
        #     label.omim_id: label.label for label in self.do_unmapped_labels
        # }
        # for omim_in_do in self.omim_in_dos:
        #     for xref in omim_in_do.xrefs:
        #         if not xref in do_by_disease_id:
        #             do_by_disease_id[xref] = omim_in_do.label
        # for label in self.do_omim_import:
        #     if label.omim_id not in do_by_disease_id:
        #         do_by_disease_id[label.omim_id] = label.label
        # Our third source for OMIM names is the CTD disease database.
        ctd_by_disease_id: Dict[str, CtdDiseaseEntry] = {}
        for ctd_disease in self.ctd_diseases:
            for disease_id in ctd_disease.alt_disease_ids:
                if disease_id not in ctd_by_disease_id:
                    ctd_by_disease_id[disease_id] = ctd_disease
        # If annotation with these names fails, Orphanet terms will be tried later.
        for record in self.mim2gene_medgen:
            disease_name = None
            disease_ids = [record.omim_id]
            if record.omim_id in hpoa_by_disease_id:
                disease_name = hpoa_by_disease_id[record.omim_id].disease_name
            elif record.omim_id in mondo_by_omim_id:
                mondo_records = mondo_by_omim_id[record.omim_id]
                disease_name = mondo_records[0].name
                disease_ids += [record.mondo_id for record in mondo_records]
            elif record.omim_id in ctd_by_disease_id:
                disease_name = ctd_by_disease_id[record.omim_id].disease_name
            elif record.omim_id in self.omim_titles:
                disease_name = self.omim_titles[record.omim_id]
            self.register_disease_assoc(
                GeneDiseaseAssociation(
                    hgnc_id=record.hgnc_id,
                    labeled_disorders=[LabeledDisorder(term_id=record.omim_id, title=disease_name)],
                    disease_ids=list(sorted(set(disease_ids))),
                    disease_name=disease_name,
                    disease_definition=None,
                    sources=[GeneDiseaseAssociationSource.OMIM],
                    confidence=ConfidenceLevel.HIGH,
                )
            )

    def handle_orpha(self):
        """Process the ORPHA data."""
        logger.info("Processing ORPHA...")
        orpha_disorders_by_orpha_id = {d.database_id: d for d in self.orpha_disorders}

        # MONDO mapping for proper aliasing.
        mondo_by_orpha_id: Dict[str, MondoDisease] = {}
        for relation in GOOD_MONDO_RELATIONS:
            for mondo_disease in self.mondo_diseases:
                for synonym in mondo_disease.synonyms:
                    if (
                        synonym.database_id.startswith("Orphanet:")
                        and relation in synonym.relation
                        and synonym.database_id not in mondo_by_orpha_id
                    ):
                        orpha_id = synonym.database_id.replace("Orphanet:", "ORPHA:")
                        mondo_by_orpha_id.setdefault(orpha_id, []).append(mondo_disease)

        for orpha_mapping in self.orpha_mappings:
            database_ids = [orpha_mapping.orpha_id]
            orpha_disorder = orpha_disorders_by_orpha_id[orpha_mapping.orpha_id]
            for mapping in orpha_disorder.mappings:
                if (
                    mapping.database_id.startswith("OMIM:")
                    and mapping.relation in GOOD_ORPHA_RELATIONS
                ):
                    database_ids.append(mapping.database_id)
            if orpha_mapping.orpha_id in mondo_by_orpha_id:
                database_ids += [
                    record.mondo_id for record in mondo_by_orpha_id[orpha_mapping.orpha_id]
                ]
            self.register_disease_assoc(
                GeneDiseaseAssociation(
                    hgnc_id=orpha_mapping.hgnc_id,
                    labeled_disorders=[
                        LabeledDisorder(
                            term_id=orpha_mapping.orpha_id, title=orpha_disorder.database_name
                        )
                    ],
                    disease_ids=list(sorted(set(database_ids))),
                    disease_name=orpha_disorder.database_name,
                    disease_definition=orpha_disorder.definition,
                    sources=[GeneDiseaseAssociationSource.ORPHANET],
                    confidence={
                        OrphaStatus.ASSESSED.value: ConfidenceLevel.HIGH,
                        OrphaStatus.NOT_YET_ASSESSED.value: ConfidenceLevel.LOW,
                    }[orpha_mapping.status],
                )
            )

    def handle_panelapp(self):  # noqa: C901
        """Process the PanelApp data."""
        RE_OMIM = r"(?:O?MIM:)?(\d\d\d\d\d\d)"
        RE_MONDO = r"MONDO:(\d+)"
        RE_ORPHA = r"ORPHA:(\d+)"

        logger.info("Processing PanelApp...")
        orpha_disorders: Dict[str, OrphaDisorder] = {
            disease.database_id: disease for disease in self.orpha_disorders
        }
        for assoc in self.panelapp_associations:
            self.panelapp_assocs.setdefault(assoc.hgnc_id, []).append(assoc)

            if assoc.confidence_level == PanelappConfidence.NONE.value:
                continue
            if assoc.phenotypes:
                for phenotype in assoc.phenotypes:
                    m_omim = re.search(RE_OMIM, phenotype, re.IGNORECASE)
                    m_mondo = re.search(RE_MONDO, phenotype, re.IGNORECASE)
                    m_orpha = re.search(RE_ORPHA, phenotype, re.IGNORECASE)
                    labeled_omim_ids: List[LabeledDisorder] = []
                    labeled_orpha_ids: List[LabeledDisorder] = []
                    omim_ids: List[LabeledDisorder] = []
                    orpha_ids: List[LabeledDisorder] = []
                    orpha_disorder: Optional[OrphaDisorder] = None
                    if m_mondo:
                        mondo_id = f"MONDO:{m_mondo.group(1)}"
                        omim_ids = list(self.mondo_to_omim.get(mondo_id, []))
                        orpha_ids = list(sorted(self.mondo_to_orpha.get(mondo_id, set())))
                        # We base the disease name and description on the first Orpha ID.
                        # If the association already exists, this one will be ignored.
                        if orpha_ids:
                            orpha_disorder = orpha_disorders.get(orpha_ids[0])
                    elif m_orpha:
                        orpha_id = f"ORPHA:{m_orpha.group(1)}"
                        orpha_ids = [orpha_id]
                        omim_ids = list(self.omim_to_orpha.get(orpha_id, []))
                        if orpha_ids:
                            orpha_disorder = orpha_disorders.get(orpha_id)
                    elif m_omim:
                        omim_id = f"OMIM:{m_omim.group(1)}"
                        if omim_id not in self.xlink_by_omim_id:  # OMIM was a gene...
                            labeled_omim_ids = omim_ids = [omim_id]
                            orpha_ids = list(sorted(self.omim_to_orpha.get(omim_id, set())))
                            # We base the disease name and description on the first Orpha ID.
                            # If the association already exists, this one will be ignored.
                            if orpha_ids:
                                orpha_disorder = orpha_disorders.get(orpha_ids[0])
                        else:
                            logger.debug(f"Unresolved phenotype {phenotype} (OMIM is gene)")
                    else:
                        logger.debug(f"Unresolved phenotype {phenotype}")

                    stripped_phenotype = (
                        phenotype.replace(f", {omim_id}", "")  # trim ", <OMIM>"
                        .replace(f" {omim_id}", "")  # trim " <OMIM>"
                        .replace(omim_id, "")  # trim "<OMIM>"
                        .replace("{", "")  # trim curly braces in "{xxx}"
                        .replace("}", "")
                    )

                    labeled_omim_ids = [
                        LabeledDisorder(
                            term_id=omim_id, title=self.omim_titles.get(omim_id, stripped_phenotype)
                        )
                        for omim_id in omim_ids
                    ]
                    labeled_orpha_ids = [
                        LabeledDisorder(
                            term_id=orpha_id, title=orpha_disorders[orpha_id].database_name
                        )
                        for orpha_id in orpha_ids
                    ]

                    if omim_ids or orpha_ids:
                        # Obtain the disease name from orpha_id/omim_id and lookup.
                        disease_name = None
                        for orpha_id in orpha_ids:
                            if orpha_id in orpha_disorders:
                                disease_name = orpha_disorders[orpha_id].database_name
                                break
                        if not disease_name:
                            for omim_id in omim_ids:
                                if omim_id in self.omim_titles:
                                    disease_name = self.omim_titles[omim_id]
                                    break
                        # Use label from PanelApp as fallback.
                        if not disease_name and m_omim:
                            disease_name = stripped_phenotype
                        # We found at least one OMIM/ORPHA ID. We will now created
                        # a new `GeneDiseaseAssociation` object from this which will
                        # be merged into any existing.
                        gene_disease_assoc = GeneDiseaseAssociation(
                            hgnc_id=assoc.hgnc_id,
                            labeled_disorders=list(
                                sorted(set(labeled_omim_ids + labeled_orpha_ids))
                            ),
                            disease_ids=list(sorted(set(omim_ids + orpha_ids))),
                            disease_name=disease_name,
                            disease_definition=orpha_disorder.definition
                            if orpha_disorder
                            else None,
                            sources=[GeneDiseaseAssociationSource.PANELAPP],
                            confidence={
                                PanelappConfidence.GREEN.value: ConfidenceLevel.HIGH,
                                PanelappConfidence.AMBER.value: ConfidenceLevel.MEDIUM,
                                PanelappConfidence.RED.value: ConfidenceLevel.LOW,
                            }[assoc.confidence_level],
                        )
                        self.register_disease_assoc(gene_disease_assoc)
            else:
                logger.debug(f"UNRESOLVED: {assoc.panel.name}")


if __name__ == "__main__":
    Integrator().run("store.pickle.gz" if DEV_MODE else None)
