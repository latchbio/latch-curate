from pathlib import Path
import uuid
import shutil
import re

import numpy as np
import anndata as ad
import scipy.sparse as sp
import docker
import textwrap
import os

from latch_curate.construct.prompts import construct_prompt, construct_instructions
from latch_curate.config import user_config
from latch_curate.constants import latch_curate_constants as lcc

def client_from_env():
    host = os.environ.get("DOCKER_HOST")
    if host is not None:
        return docker.DockerClient(base_url=host.strip())

    socket_paths = [
        "/var/run/docker.sock",  # linux
        os.path.expanduser(Path().home() / ".docker/run/docker.sock"),  # Docker Desktop on macOS newer versions
        os.path.expanduser(Path().home() / "Library/Containers/com.docker.docker/Data/docker.sock"),  # Docker Desktop on macOS older versions
    ]
    for p in socket_paths:
        if Path(p).exists():
            os.environ["DOCKER_HOST"] = f"unix://{p}"
            return docker.DockerClient(base_url=f"unix://{p}")
    raise OSError(
        "Cannot find a Docker socket. Checked:\n  " +
        "\n  ".join(socket_paths) +
        "\nEither start Docker or set DOCKER_HOST explicitly."
    )

# todo(kenny): clearly document validation criteria
# - `latch_sample_id` in obs
# - var index is ensembl
# - `gene_symbols` inv ar
# - var and obs index values are unique
# - counts are positive integers (not transformed)
# - any additiona author metadata variables are in obs prefixed with `author_`
def validate_adata(path: Path) -> None:
    adata = ad.read_h5ad(path)
    assert 'latch_sample_id' in adata.obs.columns, "missing 'latch_sample_id'"

    pattern = re.compile(r"ENS[A-Z0-9]{0,5}G\d{11}(\.\d+)?$")

    assert all(bool(pattern.match(x)) for x in adata.var_names), "var names not Ensembl IDs"

    assert 'gene_symbols' in adata.var.columns, "missing 'gene_symbols'"

    assert adata.obs_names.is_unique, "duplicate obs indices"
    assert adata.var_names.is_unique, "duplicate var indices"

    assert adata.var['gene_symbols'].is_unique, "duplicate gene symbols"

    X = adata.X
    data = X.data if sp.issparse(X) else np.asarray(X)
    assert np.all(data >= 0), "negative values in counts"
    assert np.allclose(data, np.round(data)), "non-integer counts"

    # 6) cell count roughly matches
    # expected = _extract_expected_cell_count(Path('paper_text.txt'))
    # if expected:
    #     delta = abs(adata.n_obs - expected)
    #     tol = max(100, int(expected * 0.2))
    #     assert delta <= tol, f"cell count {adata.n_obs} off from expected {expected} by > {tol}"

    for col in adata.obs.columns:
        if col.startswith('author_'):
            series = adata.obs[col]
            assert not series.isna().all(), f"{col} all NaN"
            if series.dropna().nunique() <= 1:
                raise AssertionError(f"{col} has single unique value; check source")

def construct_counts(
    data_dir: Path,
    paper_text_path: Path,
    study_metadata_path: Path,
    workdir: Path,
    #
    model: str = "o4-mini",
    max_rounds: int = 5,
) -> Path:
    workdir.mkdir(exist_ok=True)

    run_dir = workdir / uuid.uuid4().hex
    run_dir.mkdir()

    for fname in ("scrna_utils.py", "ensembl_map.json"):
        src = Path(__file__).parent / "vendor" / fname
        shutil.copy(src, run_dir / src.name)

    instructions_path = run_dir / "instructions.md"

    prompt = construct_prompt()
    instructions = construct_instructions(paper_text_path.read_text(), study_metadata_path.read_text())
    instructions_path.write_text(instructions)

    for attempt in range(1, max_rounds + 1):
        print(f"\n=== Codex attempt {attempt}/{max_rounds} ===")
        prompt_path = run_dir / f"prompt_round{attempt}.md"
        prompt_path.write_text(prompt)

        client = client_from_env()

        log_path = run_dir / f"codex_round{attempt}.log"
        with log_path.open("w") as log_f:
            codex_cmd = (
                f'codex --model {model} --approval-mode full-auto --quiet "$(cat {prompt_path.name})"'
            )
            container = client.containers.run(
                # todo(kenny): obv
                image="foobar",
                command=[
                    "/usr/bin/bash",
                    "-lc",
                    codex_cmd
                ],
                volumes={
                    str(run_dir.resolve()):  {"bind": "/workspace",      "mode": "rw"},
                    str(data_dir.resolve()): {"bind": "/workspace/data", "mode": "rw"},
                    str(instructions_path.resolve()): {"bind": "/root/.codex/instructions.md", "mode": "rw"},
                },
                environment={
                        "OPENAI_API_KEY": user_config.openai_api_key,
                        "PATH": os.environ.get("PATH", ""),
                        "PYTHONUNBUFFERED": "1",
                },
                detach=True,
                stdout=True,
                stderr=True,
                remove=True
            )

            for line in container.logs(stream=True):
                text = line.decode()
                print(text, end="")
                log_f.write(text)

        out = run_dir / 'output.h5ad'
        try:
            validate_adata(out)
            print('[ok] in-code validation passed')
            break
        except Exception as err:
            print(f'[retry] validation failed: {err}')
            excerpt = textwrap.indent(str(err), '> ')
            prompt += f"\n\n# Validation failure: {excerpt}\n"
    else:
        raise RuntimeError("Maximum rounds reached without passing tests.")

    artefact_src = run_dir / "output.h5ad"
    if not artefact_src.exists():
        raise RuntimeError("output.h5ad missing despite test success.")

    artefact_dest = workdir / lcc.construct_counts_adata_name
    shutil.move(artefact_src, artefact_dest)
    print(f"anndata at {artefact_dest}")
    return artefact_dest
