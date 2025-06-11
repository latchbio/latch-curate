import json
from textwrap import dedent
from pathlib import Path

def _quote_sql(val: str | int | None) -> str:
    if val is None:
        return "NULL"
    if isinstance(val, (int, float)):
        return str(val)
    escaped = str(val).replace("'", "''")
    return f"'{escaped}'"

def build_sql(data: dict, *, ldata_node_id: int, provider: str, price: int, agreement_id: int) -> str:

    info = data["info"]
    info_row = [
        provider,
        price,
        ldata_node_id,
        info["display_name"].strip(),
        info["description"].strip(),
        info["paper_title"].strip(),
        info["paper_url"].strip(),
        info["data_url"].strip().format(gse_id=info["data_external_id"].strip()),
        agreement_id,
        info["cell_count"],
        info["data_external_id"].strip(),
    ]
    info_sql = dedent(f"""
        INSERT INTO app_public.dataset_info
            (provider, price, ldata_node_id, display_name, description,
             paper_title, paper_url, data_url, agreement_id, cell_count,
             data_external_id)
        VALUES
            ({", ".join(_quote_sql(v) for v in info_row)});
    """).strip()

    tag_rows = []
    for t in data["tags"]:
        tag_rows.append(
            f"({ldata_node_id}, "
            f"{_quote_sql(t['metadata_type'])}, "
            f"{_quote_sql(t['value'])}, "
            f"{_quote_sql(t.get('ontology_id') or '')})"
        )
    tag_values_block = ",\n".join(tag_rows)

    tags_sql = dedent(f"""
        INSERT INTO app_public.dataset_info_tag (dataset_id, key, value, ontology_id)
        SELECT di.id,
               t.key::app_public.dataset_metadata_type,
               t.value,
               t.ontology_id
        FROM (VALUES
        {tag_values_block}
        ) AS t(ldata_node_id, key, value, ontology_id)
        JOIN app_public.dataset_info AS di
          ON di.ldata_node_id = t.ldata_node_id;
    """).strip()

    return info_sql + "\n\n" + tags_sql + "\n"

def build_publish_queries(ldata_node_id: int, build_data_file: Path, build_queries_file: Path) -> None:

    data = json.loads(build_data_file.read_text())

    provider = "orion"
    price = data["info"]["cell_count"] * 0.01

    sql = build_sql(
        data,
        ldata_node_id=ldata_node_id,
        provider=provider,
        price=price,
        agreement_id=1,
    )

    build_queries_file.write_text(sql)
    print(f"Wrote {build_queries_file}")
