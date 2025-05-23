import io
import base64
import matplotlib.pyplot as plt
import pandas as pd

import openai

from latch_curate.config import OPENAI_KEY

def _fig_to_base64(obj) -> str:

    if isinstance(obj, (list, tuple)) and obj:
        obj = obj[0]

    if hasattr(obj, "savefig"):
        fig = obj
    elif hasattr(obj, "figure"):
        fig = obj.figure
    elif hasattr(obj, "fig"):
        fig = obj.fig
    else:
        raise TypeError(f"Cannot extract Figure from object of type {type(obj)}")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")

def print_full_df_string(df: pd.DataFrame) -> str:
    return df.to_string(max_rows=None, max_cols=None, line_width=None)


def prompt_model(prompt: str, model: str = "o4-mini") -> str:
    openai.api_key = OPENAI_KEY

    response = openai.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
    )

    return response.choices[0].message.content.strip()
