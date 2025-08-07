from pathlib import Path
from dataclasses import dataclass

import pandas as pd
from anndata import AnnData

from latch_curate.harmonize.schema import parse_metadata_yaml, Vocab
from latch_curate.publish.types import Tag
from latch_curate.lint.vocab.cl import validate_cl
from latch_curate.lint.vocab.mondo import validate_mondo
from latch_curate.lint.vocab.uberon import validate_uberon
from latch_curate.lint.vocab.sample_site import validate_sample_site
from latch_curate.lint.vocab.sequencing_platform import validate_efo
from latch_curate.lint.vocab.common import unknown_val


@dataclass
class ValidationError:
    column: str
    value: str
    error_type: str
    message: str


@dataclass
class ValidationResult:
    is_valid: bool
    errors: list[ValidationError]
    tags: list[Tag]
    validated_columns: set[str]


class MetadataValidator:

    ontology_validators = {
        "cl": validate_cl,
        "mondo": validate_mondo, 
        "uberon": validate_uberon,
        "sample_site": validate_sample_site,
        "efo": validate_efo,
    }


    vocab_to_tag_type = {
        "cl": "cell_type",
        "mondo": "disease", 
        "uberon": "tissue",
        "efo": "assay",
        "sample_site": "sample_site",
    }

    def __init__(self, schema_path: Path):
        self.variable_definitions = parse_metadata_yaml(schema_path)
        self.var_def_map = {var.name: var for var in self.variable_definitions}

    def validate_obs_metadata(self, adata: AnnData) -> ValidationResult:
        errors = []
        tags = []
        validated_columns = set()

        for var_def in self.variable_definitions:
            column_name = var_def.name

            if column_name not in adata.obs.columns:
                errors.append(ValidationError(
                    column=column_name,
                    value="",
                    error_type="missing_column",
                    message=f"Required column '{column_name}' is missing from obs"
                ))
                continue

            validated_columns.add(column_name)
            unique_values = set(adata.obs[column_name].dropna().astype(str))

            column_errors, column_tags = self._validate_column_values(
                column_name, unique_values, var_def.vocab
            )
            errors.extend(column_errors)
            tags.extend(column_tags)

        is_valid = len(errors) == 0

        return ValidationResult(
            is_valid=is_valid,
            errors=errors,
            tags=tags,
            validated_columns=validated_columns
        )

    def _validate_column_values(
        self, column_name: str, values: set[str], vocab: Vocab
    ) -> tuple[list[ValidationError], list[Tag]]:
        errors = []
        tags = []

        for value in values:
            if value == unknown_val or pd.isna(value) or value == "":
                continue

            error, tag = self._validate_single_value(column_name, value, vocab)
            if error:
                errors.append(error)
            if tag:
                tags.append(tag)

        return errors, tags

    def _validate_single_value(
        self, column_name: str, value: str, vocab: Vocab
    ) -> tuple[ValidationError, Tag]:

        if vocab.type == "ontology":
            return self._validate_ontology_value(column_name, value, vocab)
        elif vocab.type == "custom":
            return self._validate_custom_value(column_name, value, vocab)
        elif vocab.type == "uncontrolled":
            return self._validate_uncontrolled(column_name, value, vocab)
        else:
            return ValidationError(
                column=column_name,
                value=value,
                error_type="unknown_vocab_type",
                message=f"Unknown vocabulary type: {vocab.type}"
            ), None

    def _validate_ontology_value(
        self, column_name: str, value: str, vocab: Vocab
    ) -> tuple[ValidationError, Tag]:
        if not vocab.name:
            return ValidationError(
                column=column_name,
                value=value,
                error_type="missing_ontology_name",
                message="Ontology vocabulary missing 'name' field"
            ), None

        validator = self.ontology_validators.get(vocab.name)
        if not validator:
            return ValidationError(
                column=column_name,
                value=value,
                error_type="unsupported_ontology",
                message=f"Unsupported ontology: {vocab.name}"
            ), None

        is_valid = validator(value)
        if not is_valid:
            return ValidationError(
                column=column_name,
                value=value,
                error_type="invalid_ontology_term",
                message=f"'{value}' is not a valid {vocab.name} term"
            ), None


        tag_type = self.vocab_to_tag_type.get(vocab.name)
        if tag_type and "/" in value:
            name, obo_id = value.split("/", 1)
            return None, Tag(metadata_type=tag_type, name=name, obo_id=obo_id)

        return None, None

    def _validate_custom_value(
        self, column_name: str, value: str, vocab: Vocab
    ) -> tuple[ValidationError, Tag]:
        if not vocab.values:
            return ValidationError(
                column=column_name,
                value=value,
                error_type="missing_custom_values",
                message="Custom vocabulary missing 'values' field"
            ), None

        if value not in vocab.values:
            return ValidationError(
                column=column_name,
                value=value,
                error_type="invalid_custom_value",
                message=f"'{value}' not in allowed values: {vocab.values}"
            ), None

        return None, None

    def _validate_uncontrolled(
        self, column_name: str, value: str, vocab: Vocab
    ) -> tuple[ValidationError, Tag]:
        if not value or value.strip() == "":
            return ValidationError(
                column=column_name,
                value=value,
                error_type="empty_uncontrolled",
                message="Uncontrolled field cannot be empty"
            ), None

        return None, None

    def print_validation_report(self, result: ValidationResult):
        print("=" * 60)
        print("METADATA VALIDATION REPORT")
        print("=" * 60)

        print(f"\nValidation Status: {'✓ PASSED' if result.is_valid else '✗ FAILED'}")
        print(f"Columns Validated: {len(result.validated_columns)}")
        print(f"Errors Found: {len(result.errors)}")
        print(f"Tags Extracted: {len(result.tags)}")

        if result.validated_columns:
            print("\nValidated Columns:")
            for col in sorted(result.validated_columns):
                print(f"  - {col}")

        if result.errors:
            print("\nValidation Errors:")
            for error in result.errors:
                print(f"  ✗ [{error.column}] {error.message}")
                if error.value:
                    print(f"    Value: '{error.value}'")

        if result.tags:
            print("\nExtracted Tags:")
            for tag in result.tags:
                print(f"  - {tag.metadata_type}: {tag.name} ({tag.obo_id})")

        print("=" * 60)


def validate_harmonized_metadata(
    adata: AnnData, schema_path: Path
) -> ValidationResult:
    validator = MetadataValidator(schema_path)
    result = validator.validate_obs_metadata(adata)
    validator.print_validation_report(result)
    return result
