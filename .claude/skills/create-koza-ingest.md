---
name: create-koza-ingest
description: Add a new Koza 2.x ingest to this repository
version: 2.0.0
triggers:
  - "create a koza ingest"
  - "new koza ingest"
  - "add ingest"
---

# Create Koza Ingest Skill

You are helping to add a new Koza 2.x ingest to an existing multi-ingest repository following Monarch Initiative patterns.

## Workflow Overview

Follow these steps in order:

### 1. Gather Ingest Details

**Ask the user for:**
- Ingest name (snake_case, e.g., "gene_to_disease")
- Data source URL (direct download link)
- File format (TSV, CSV, JSON, etc.)
- Brief description of what the ingest produces

### 2. Update download.yaml with New Entry

Add a new entry to `download.yaml` (at repo root):

```yaml
downloads:
  # Existing entries above...
  - url: "https://example.org/data/file.tsv"
    local_name: "data/<ingest_name>.tsv"
```

**IMPORTANT:** Append to the existing downloads list, do not replace it.

### 3. Download the Data

```bash
just download
```

This downloads all files in download.yaml to the `data/` directory.

### 4. Inspect the Downloaded Data

Examine the data file structure:
- Check number of rows: `wc -l data/<filename>`
- View first rows: `head -20 data/<filename>`
- Look for:
  - Header rows (and whether they have comment characters)
  - Column delimiters (tab, comma, etc.)
  - Comment lines
  - Data patterns

### 5. Create <ingest_name>.yaml Transform Config

Create `src/<ingest_name>.yaml`:

```yaml
name: "<ingest_name>"
description: "Description of what this ingest produces"
metadata:
  ingest_title: "<Human Readable Title>"
  ingest_url: "https://source-url.org"
reader:
  format: "csv"
  files:
    - "../../data/<filename>"
  delimiter: "\t"
transform_code: "<ingest_name>.py"
transform_mode: "flat"
writer:
  min_node_count: 0
  min_edge_count: 0
```

#### Handling Comment Characters with Headers

When a file has comment lines AND a header row with a comment character prefix:

```yaml
reader:
  format: "csv"
  files:
    - "../../data/<filename>"
  delimiter: "\t"
  comment_char: "#"  # Skip lines starting with #
  header_mode: 1  # After filtering comments, line 1 is the header (1-indexed)
  header_prefix: "#"  # Strip "#" prefix from header line column names
```

**Example file:**
```
# See end of file for documentation
# Phenotype	Gene Symbols	MIM Number	Cyto Location
17,20-lyase deficiency, isolated, 202110 (3)	CYP17A1, CYP17, P450C17	609300	10q24.32
```

With this config, you get clean column names: `Phenotype`, `Gene Symbols`, `MIM Number`, `Cyto Location`

**Pattern:** `comment_char` + `header_mode` + `header_prefix` is a powerful combination for messy data files.

### 6. Test Configuration by Printing Records

Create a temporary transform to verify the config is correct.

Create `src/<ingest_name>.py`:

```python
from typing import Any

import koza
from koza import KozaTransform
from biolink_model.datamodel.pydanticmodel_v2 import Entity, Association


@koza.transform_record()
def transform_record(koza_transform: KozaTransform, row: dict[str, Any]) -> list[Entity | Association]:
    """Test transform - prints records to verify config."""
    print(f"Row keys: {list(row.keys())}")
    print(f"Row data: {row}")
    print("-" * 80)
    return []  # Return empty list for now
```

Run with row limit:

```bash
uv run koza transform --source src/<ingest_name>.yaml --row-limit 3
```

Verify that:
- Column names are correct
- Data is parsed properly
- No unexpected characters or issues

### 7. Research Biolink Model for Appropriate Predicates and Association Types

**CRITICAL:** Cross-reference data documentation with Biolink model.

#### Load Biolink Model Using SchemaView

**The pattern for researching Biolink model predicates and classes:**

```python
from importlib.resources import files
from linkml_runtime.utils.schemaview import SchemaView

# Load Biolink model from installed biolink_model package
biolink_yaml = files('biolink_model.schema') / 'biolink_model.yaml'
sv = SchemaView(str(biolink_yaml))

# Look up predicates (slots in LinkML terminology)
slot = sv.get_slot("causes")
print(f"Description: {slot.description}")
print(f"Parent: {slot.is_a}")
print(f"Exact mappings: {slot.exact_mappings}")  # RO terms here!
print(f"Domain: {slot.domain}")
print(f"Range: {slot.range}")

# Look up association classes
cls = sv.get_class("causal gene to disease association")
print(f"Description: {cls.description}")
print(f"Parent: {cls.is_a}")

# List all available predicates/slots
all_slots = sv.all_slots()
print([s for s in all_slots if 'predispose' in s.lower()])
print([s for s in all_slots if 'cause' in s.lower()])
```

**Why use SchemaView instead of docs?**
- Get exact mappings to RO terms
- See parent/child relationships between predicates
- Discover related predicates you might not know about
- Access programmatically for validation scripts
- Always reflects the installed biolink model version

**Map data semantics to Biolink:**
1. Check predicate **descriptions** - does it match your relationship?
2. Check predicate **exact_mappings** - what RO term does it map to?
3. Check predicate **parent** (is_a) - understand hierarchy
4. Determine which **association class** to use
5. Match data confidence levels/markers to appropriate predicates

### 8. Validate Predicate Choices with Relation Ontology (RO) Terms

**CRITICAL:** When choosing between similar predicates, research the actual RO terms to justify your choice.

#### Check Biolink Mappings to RO

Use SchemaView to see what RO terms predicates map to:

```python
from importlib.resources import files
from linkml_runtime.utils.schemaview import SchemaView

biolink_yaml = files('biolink_model.schema') / 'biolink_model.yaml'
sv = SchemaView(str(biolink_yaml))

# Compare candidate predicates
predisposes = sv.get_slot("predisposes to condition")
print(f"Description: {predisposes.description}")
print(f"Exact mappings: {predisposes.exact_mappings}")
print(f"Broad mappings: {predisposes.broad_mappings}")

contributes = sv.get_slot("contributes to")
print(f"Description: {contributes.description}")
print(f"Exact mappings: {contributes.exact_mappings}")
print(f"Narrow mappings: {contributes.narrow_mappings}")
```

#### Search Ontology Lookup Service (OLS) for RO Terms

Use the OLS API (available at `https://www.ebi.ac.uk/ols4/api/mcp`) to look up RO terms and their definitions.

**Use WebFetch to search OLS:**

```bash
# Search for concept-related RO terms
WebFetch: https://www.ebi.ac.uk/ols4/api/search?q=susceptibility&ontology=ro
Prompt: "List all matching RO terms with their IDs, labels, and definitions"

WebFetch: https://www.ebi.ac.uk/ols4/api/search?q=predisposition&ontology=ro
Prompt: "List all matching RO terms with their IDs, labels, and definitions"

# Look up specific RO term by ID
WebFetch: https://www.ebi.ac.uk/ols4/api/search?q=RO:0019501
Prompt: "Get full details for this RO term including definition, synonyms, and related terms"

WebFetch: https://www.ebi.ac.uk/ols4/api/search?q=RO:0002326
Prompt: "Get full details for this RO term including definition, synonyms, and related terms"
```

**Example findings from OLS:**
- **RO:0019501** "confers susceptibility to condition": "Relates a gene to condition, such that a variation in this gene predisposes to the development of a condition"
- **RO:0002326** "contributes to": General contribution relationship (e.g., enzyme subunits contributing to enzyme activity)
- **RO:0003303** "causes": Direct causation
- **RO:0004015** "is causal susceptibility factor for": Necessary but not sufficient for disease development

**Key principle:** If RO has a **specific term** for your relationship type (e.g., RO:0019501 for susceptibility), that indicates the relationship is semantically distinct and should use a specific predicate rather than a generic one.

**Pattern:** Strong predicate justification = Biolink description + RO term definition + data examples

### 9. Create Unit Tests with Actual Data (BEFORE Implementing Transform)

**CRITICAL: Tests come BEFORE implementation!**

Create test file `tests/test_<ingest_name>.py` with fixtures using actual data rows:

```python
import pytest
from biolink_model.datamodel.pydanticmodel_v2 import GeneToDiseaseAssociation
from <ingest_name> import transform_record


@pytest.fixture
def standard_row_entities():
    """Standard row with expected behavior."""
    row = {
        "Column1": "value1",
        "Column2": "value2",
    }
    # Pass None for koza_transform since we don't use maps/lookups
    # (Note: pass actual koza_transform if transform uses get_map())
    return transform_record(None, row)


def test_standard_case(standard_row_entities):
    """Test standard case."""
    assert standard_row_entities
    assert len(standard_row_entities) == 1

    association = standard_row_entities[0]
    assert isinstance(association, GeneToDiseaseAssociation)
    assert association.subject == "EXPECTED:ID"
    assert association.object == "EXPECTED:ID"
    assert association.predicate == "biolink:expected_predicate"
```

**Testing pattern notes:**
- Call `transform_record(None, row)` directly for transforms without maps
- Pass actual `koza_transform` object if transform uses `get_map()` for lookups
- Create fixtures for ALL edge cases found in documentation research
- Tests should FAIL initially (since transform not implemented yet)

### 10. Implement <ingest_name>.py

Now implement the actual transform logic based on:
- Test expectations
- Data format documentation
- Biolink model mappings

```python
import uuid
import re
from typing import Any

import koza
from koza import KozaTransform
from biolink_model.datamodel.pydanticmodel_v2 import (
    CausalGeneToDiseaseAssociation,
    CorrelatedGeneToDiseaseAssociation
)


@koza.transform_record()
def transform_record(koza_transform: KozaTransform, row: dict[str, Any]) -> list:
    """Transform data row into Biolink entities/associations."""

    # Parse and extract IDs
    # Handle special cases (brackets, question marks, etc.)
    # Determine appropriate association type and predicate
    # Create and return associations

    return [association]
```

### 11. Update TRANSFORMS in justfile

Add the new ingest to the TRANSFORMS variable in `justfile`:

```just
# Find this line and add your new ingest:
TRANSFORMS := "existing_ingest another_ingest <ingest_name>"
```

This allows `just transform` and `just test` to include the new ingest.

### 12. Run Tests and Verify

```bash
# Run tests for all ingests
just test

# Or run tests for just your new ingest
uv run pytest tests/test_<ingest_name>.py -v

# Run full transform for your ingest
uv run koza transform --source src/<ingest_name>.yaml

# Check output
ls -lh output/
```

### 13. Update min_edge_count/min_node_count

Based on actual output, update the YAML config:

```yaml
writer:
  min_node_count: 0  # Set based on expected output
  min_edge_count: <number>  # Set slightly less than actual count for buffer
```

## Key Patterns Summary

1. **Always configure YAML first, test with printed records before writing transform code**
2. **Research data format documentation thoroughly**
3. **Cross-reference with Biolink model using SchemaView**
4. **Validate predicate choices with RO terms from OLS** - if RO has a specific term for your relationship, use a specific predicate
5. **Write tests BEFORE implementation, using actual data rows**
6. **Test pattern: `transform_record(None, row)` unless using maps, then pass `koza_transform`**
7. **Comment character + header_mode + header_prefix pattern for messy headers**
8. **Update TRANSFORMS in justfile when adding new ingests**
9. **Use `just` commands: `just download`, `just test`, `just transform`**

## Files Created/Modified for Each New Ingest

| File | Action |
|------|--------|
| `download.yaml` | Append new download entry |
| `src/<ingest_name>.yaml` | Create new transform config |
| `src/<ingest_name>.py` | Create new transform code |
| `tests/test_<ingest_name>.py` | Create new test file |
| `justfile` | Update TRANSFORMS variable |
