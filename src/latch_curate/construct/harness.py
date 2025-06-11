from dataclasses import dataclass
from pathlib import Path
import shutil

import docker
import os
import anndata as ad
from anndata import AnnData
import pandas as pd
import scipy.sparse as sp
from textwrap import dedent, indent, shorten
import json

from latch_curate.llm_utils import prompt_model
from latch_curate.constants import latch_curate_constants as lcc
from latch_curate.construct.validate import validate_counts_object
from latch_curate.utils import _df_to_html, write_html_report
from latch_curate.construct.prompts import build_construct_counts_prompt, build_construct_counts_instructions, build_get_target_cell_count_prompt
from latch_curate.config import user_config

OPENAI_SYSTEM_PROMPT_LEN = 2512
OPENAI_MAX_PROMPT_LEN = 1048576

def get_system_memory() -> int:
    pages = os.sysconf('SC_PHYS_PAGES')
    page_size = os.sysconf('SC_PAGE_SIZE')
    total_bytes = pages * page_size
    return total_bytes

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

def build_construct_report_html(
        adata: AnnData,
        checks: list[tuple[str, str]],
) -> str:
    title = "Construct-Counts Report"

    check_df = pd.DataFrame(checks, columns=["status", "check"])
    check_df["status"] = check_df["status"].map(
        {"pass": "✅", "fail": "❌"}
    )
    validation_html = _df_to_html(check_df)

    nnz = adata.X.nnz if sp.issparse(adata.X) else (adata.X != 0).sum()
    stats_df = pd.DataFrame(
        {
            "metric": [
                "# cells",
                "# genes",
                "non-zero counts",
                "sparsity %"],
            "value": [
                adata.n_obs,
                adata.n_vars,
                nnz,
                100 * (1 - nnz / (adata.n_obs * adata.n_vars))
            ],
        }
    )
    stats_html = _df_to_html(stats_df)

    obs_head_html = _df_to_html(adata.obs.head())

    body = dedent(
        f"""
        <h2>Validation suite</h2>
        {validation_html}

        <h2>AnnData statistics</h2>
        {stats_html}

        <h2>obs (head)</h2>
        {obs_head_html}
        """
    )

    return dedent(
        f"""<!DOCTYPE html>
        <html><head><meta charset='utf-8'>
          <title>{title}</title>
          <style>
            body {{ font-family:sans-serif; max-width:900px; margin:auto; }}
            h2 {{ border-bottom:1px solid #ccc; margin-top:2em; }}
            table {{ border-collapse:collapse; width:100%; margin-bottom:1em; }}
            th,td {{ border:1px solid #ddd; padding:4px 8px; text-align:center; }}
            .quant {{ font-size:0.8em; }}
          </style>
        </head><body>
          <h1>{title}</h1>
          {body}
        </body></html>"""
    )


@dataclass
class AgentStep:
    id: str
    kind: str
    subtype: str
    call_id: str | None
    # todo(kenny): structure metadata
    meta: dict[str, any]

def _short(s: any, trunc: int, ellip: str) -> str:
    if s is None:
        return ""
    if not isinstance(s, str):
        s = json.dumps(s, ensure_ascii=False)
    s = s.replace("\n", " ")
    return shorten(s, trunc, placeholder=ellip)


def pretty_step(step: AgentStep, width: int = 500) -> str:
    TRUNC = width - 40
    ellip = "…"

    prefix = step.id[:8]

    if step.kind == "reasoning":
        duration_ms = step.meta.get("duration_ms", 0)
        if duration_ms is not None:
            dur = step.meta.get("duration_ms", 0) / 1000
            summary_data = step.meta.get('summary')
            if type(summary_data) is not list:
                summary_data = summary_data.get('text')
            return f"{prefix} THINK  {dur:6.2f}s {summary_data}"
        return f"{prefix} THINK"

    if step.kind == "function_call": 
        args = _short(step.meta.get("arguments"), TRUNC, ellip)
        return f"{prefix} CALL   {step.subtype:<15} {args}"

    if step.kind == "function_call_output": 
        try:
            data = json.loads(step.meta.get('output'))
            output = data['output']
            metadata = data['metadata']
            exit_code = metadata['exit_code']
            if exit_code == 137:
                # todo(kenny): don't raise in here
                raise SystemError("Agent container OOM. Increase resources.")
            dur = metadata['duration_seconds']
        except Exception:
            return
        return f"{prefix} OUTPUT {dur:6.2f}s code {exit_code} \n\t>>>\t{_short(output, TRUNC, ellip)}"

    if step.kind == "message":
        if prefix == "<no-id>":
            return

        content = step.meta.get("content")
        txt = (
            " ".join(c.get("text", "") for c in content)
            if isinstance(content, list) else str(content)
        ) 
        return f"{prefix} MSG    {txt}"

    return f"{prefix} {step.kind.upper():7} {step.subtype}"

def parse_codex_line(raw: str) -> AgentStep | None:

    raw = raw.strip()
    if not raw:
        return None
    try:
        rec = json.loads(raw)
    except json.JSONDecodeError:
        return None

    kind = rec.get("type") or "<unknown>"
    call_id = rec.get("call_id")
    subtype, meta = kind, {}

    if kind == "reasoning":
        meta = {"duration_ms": rec.get("duration_ms"),
                "summary": rec.get("summary")}
    elif kind == "function_call":
        subtype = rec.get("name") or "<unknown-fn>"
        meta = {"status": rec.get("status"),
                "arguments": rec.get("arguments")}
        call_id = rec.get("call_id")
    elif kind == "function_call_output":
        subtype = "output"
        meta = {"output": rec.get("output")}
    elif kind == "message":
        subtype = "assistant_message"
        meta = {"content": rec.get("content")}

    return AgentStep(id=rec.get("id") or rec.get("call_id") or "<no-id>", kind=kind, subtype=subtype, call_id=call_id, meta=meta)


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

    for fname in ("scrna_utils.py", "ensembl_map.json"):
        src = Path(__file__).parent / "vendor" / fname
        shutil.copy(src, workdir / src.name)

    instructions_path = workdir / "instructions.md"
    construct_counts_prompt_path = workdir / "construct_counts_prompt.md"

    paper_text = paper_text_path.read_text()
    study_metadata = study_metadata_path.read_text()

    get_target_cell_count_prompt = build_get_target_cell_count_prompt(paper_text, study_metadata)

    while True:
        print("Requesting target cell count from language model")
        message_resp_json, _ = prompt_model([{"role": "user", "content": get_target_cell_count_prompt}])
        try:
            data = json.loads(message_resp_json)
            target_cell_count = int(data["target_cell_count"])
            reasoning = data["reasoning"]
            print(f" target cell count >> {target_cell_count}")
            print(f" reasoning >> {reasoning}")
            break
        except Exception:
            print(f">>> resp: {message_resp_json}")
            print("Malformed JSON output model. Trying again.")
            continue

    construct_counts_prompt = build_construct_counts_prompt(target_cell_count)
    instructions = build_construct_counts_instructions(paper_text, study_metadata)

    total_size = len(instructions) + len(construct_counts_prompt) + OPENAI_SYSTEM_PROMPT_LEN
    if len(instructions) + len(construct_counts_prompt) + OPENAI_SYSTEM_PROMPT_LEN > OPENAI_MAX_PROMPT_LEN:
        raise ValueError(
            f"Input too large ({total_size}/{OPENAI_MAX_PROMPT_LEN} chars). "
            "Please shorten your paper or supplementary files in the `download`"
            " folder by removing unnecessary text, like references, links, etc."
        )
    instructions_path.write_text(instructions)
    construct_counts_prompt_path.write_text(construct_counts_prompt)

    for attempt in range(1, max_rounds + 1):
        print(f"\n=== Agential count matrix construction attempt {attempt}/{max_rounds} ===")

        client = client_from_env()
        system_memory_bytes = get_system_memory()
        docker_mem_limit = int(system_memory_bytes * 0.8) // (1024 ** 3)
        print(f"{system_memory_bytes} bytes of system memory. Setting agent limit to {docker_mem_limit} GB.")
        log_path = workdir / f"codex_round{attempt}.log"
        with log_path.open("w") as log_f:
            codex_cmd = (
                f'codex --model {model} --approval-mode full-auto --quiet "$(cat {construct_counts_prompt_path.name})"'
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
                    str(workdir.resolve()):  {"bind": "/workspace",      "mode": "rw"},
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
                remove=True,
                mem_limit=f"{docker_mem_limit}g"
            )

            steps = []
            for line in container.logs(stream=True):
                text = line.decode()
                step = parse_codex_line(text)
                if step is not None:
                    steps.append(step)
                    print(pretty_step(step))
                log_f.write(text)

        out = workdir / 'output.h5ad'
        try:
            adata = ad.read_h5ad(out)
            validation_log = validate_counts_object(adata)
            print('[ok] in-code validation passed')
            break
        except Exception as err:
            print(f'[retry] validation failed: {err}')
            excerpt = indent(str(err), '> ')
            construct_counts_prompt += f"\n\n# Validation failure: {excerpt}\n"
    else:
        html = build_construct_report_html(adata, validation_log)
        write_html_report(html, workdir, lcc.construct_counts_report_name)
        raise RuntimeError("Maximum rounds reached without passing tests.")

    artefact_src = workdir / "output.h5ad"
    if not artefact_src.exists():
        raise RuntimeError("output.h5ad missing despite test success.")

    html = build_construct_report_html(adata, validation_log)
    write_html_report(html, workdir, lcc.construct_counts_report_name)

    artefact_dest = workdir / lcc.construct_counts_adata_name
    shutil.move(artefact_src, artefact_dest)
    print(f"anndata at {artefact_dest}")
    return artefact_dest
