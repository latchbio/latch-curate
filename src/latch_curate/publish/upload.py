from pathlib import Path
import yaml
import requests
import click
import pandas as pd

from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.config import user_config
from latch_cli.services.cp.main import cp as latch_cp
from latch.ldata._transfer.node import get_node_data
from latch.ldata._transfer.progress import Progress


def upload_dataset(
    latch_dest: str | None = None,
    curator_id: int | None = None,
    version: str | None = None,
    curator_dataset_id: str | None = None
) -> None:

    publish_workdir = Path.cwd() / lcc.publish_workdir_name
    build_file = publish_workdir / lcc.publish_build_info_file_name

    if not build_file.exists():
        print(f"Error: Build file not found at {build_file}")
        print("Run 'latch-curate publish build' first to generate build data.")
        return

    with open(build_file, 'r') as f:
        publish_data = yaml.safe_load(f)

    print("Dataset Upload")
    print(f"Paper Title: {publish_data['info']['paper_title']}")
    print(f"Cell Count: {publish_data['info']['cell_count']:,}")
    print(f"Validation Status: {publish_data['validation']['status']}")
    print(f"Tags: {len(publish_data['tags'])} extracted")

    if latch_dest is None:
        latch_dest = click.prompt("Enter destination path in Latch Data")

    if curator_id is None:
        curator_id = click.prompt("Enter curator organization ID", type=int)

    if version is None:
        default_version = "v0.1.0"
        version = click.prompt("Enter dataset version", default=default_version)

    if curator_dataset_id is None:
        default_id = publish_data['info']['data_external_id']
        curator_dataset_id = click.prompt("Enter curator dataset ID", default=default_id)

    if not version.startswith('v'):
        version = f'v{version}'

    publish_data['curator'] = {
        'curator_id': curator_id,
        'version': version,
        'curator_dataset_id': curator_dataset_id,
        'upload_timestamp': pd.Timestamp.now().isoformat()
    }

    with open(build_file, 'w') as f:
        yaml.safe_dump(publish_data, f, default_flow_style=False, indent=2)

    print("Uploading dataset...")
    print(f"Curator ID: {curator_id}")
    print(f"Version: {version}")
    print(f"Dataset ID: {curator_dataset_id}")

    full_latch_dest = f'{latch_dest}/{curator_dataset_id}'

    latch_cp([str(publish_workdir.resolve())], full_latch_dest, progress=Progress.tasks, verbose=False, expand_globs=False)
    res = get_node_data(full_latch_dest)
    node_id = res.data[full_latch_dest].id
    print(f'Retrieved node ID {node_id} for {full_latch_dest}')

    publish_data['curator']['ldata_node_id'] = node_id
    with open(build_file, 'w') as f:
        yaml.safe_dump(publish_data, f, default_flow_style=False, indent=2)

    print(f'Uploading dataset to active workspace {node_id}')

    payload = {
        "ldata_node_id": int(node_id),
        "curator_id": curator_id,
        "curator_dataset_id": curator_dataset_id,
        "version": version,
        "paper_title": publish_data['info']['paper_title'],
        "paper_url": publish_data['info']['paper_url'],
        "data_url": publish_data['info']['data_url'],
        "description": publish_data['info']['description'],
        "tags": [
            {
                "key": tag['metadata_type'],
                "value": tag['value'],
                "ontology_id": tag['ontology_id']
            }
            for tag in publish_data['tags']
        ],
        "session_id": -1,
        "workspace_id": int(user_config.workspace_id)
    }

    try:
        resp = requests.post(
            f"{lcc.nucleus_url}/publish-dataset",
            json=payload,
            headers={"Authorization": f"Latch-SDK-Token {user_config.token}"},
            timeout=60
        )
        resp.raise_for_status()

        result = resp.json()['data']
        if result['success']:
            print("Upload complete!")
            print(f"Family ID: {result['dataset_family_id']}")
            print(f"Dataset ID: {result['dataset_info_id']}")
        else:
            print(f"Upload failed: {result['message']}")

    except Exception as e:
        print(f"Error calling publish endpoint: {e}")
        print("Dataset uploaded to Latch Data but not registered in database")
